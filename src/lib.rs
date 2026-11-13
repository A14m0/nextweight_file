//! A speedier weight file library.
//!
//! A weight file can be loaded whole with [`NextWeightFile`], or walked a
//! polyid at a time with [`NwtReader`] and [`NwtWriter`] when it is too large
//! to sit in memory comfortably. [`NextWeightFile::convert_weight_file`]
//! converts a NetCDF weight file straight to disk the same way, holding no
//! more than one polyid's worth of grid data at a time.
//!
//! The on-disk layout is:
//!
//! ```text
//! "NEWT" | u64 json_len | u64 num_polyids | u64 lat_len | u64 lon_len
//!        | u64 json_offset | u64 lookup_offset
//! JSON   | the attributes and polyid names
//! lookup | (u64 data offset, u64 point count) per polyid
//! data   | (u32 lat_idx, u32 lon_idx, f32 lat, f32 lon, f32 weight) per point
//! ```

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::mem::size_of;
use std::path::{Path, PathBuf};

use netcdf::AttributeValue;

/// Magic bytes every weight file starts with
const MAGIC: &[u8; 4] = b"NEWT";
/// The fixed header: the magic bytes followed by six `u64` fields
const HEADER_LEN: usize = 4 + 6 * size_of::<u64>();
/// One lookup table entry: where a polyid's points start and how many it has
const LOOKUP_ENTRY_LEN: usize = 2 * size_of::<u64>();
/// One gridpoint on disk: two indices, two coordinates and a weight
const RECORD_LEN: usize = 5 * size_of::<u32>();
/// How much of a file to buffer while streaming through it
const STREAM_BUFFER_LEN: usize = 1 << 20;

/// A single gridpoint: latitude index, longitude index, latitude, longitude
/// and the weight itself
pub type GridPoint = (u32, u32, f32, f32, f32);

#[derive(Debug)]
pub struct NextWeightFile {
    json_data: JsonData,
    lat_len: u64,
    lon_len: u64,
    polyid_gridpoints: Vec<PolyidEntry>,
    lookup_table: Vec<(u64, u64)>,
}

#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct JsonData {
    global_attrs: Vec<(String, String)>,
    // a BTreeMap (not a HashMap) so the serialized header, and with it the
    // whole .nwt file, is byte-for-byte reproducible between runs
    per_variable_attrs: BTreeMap<String, Vec<(String, String)>>,
    polyids: Vec<String>
}

#[derive(Debug, Clone)]
#[repr(C)]
pub struct PolyidEntry {
    // lat_idx, lon_idx, lat_val, lon_val, data_value
    pub data: Vec<GridPoint>
}

impl NextWeightFile {
    /// opens a NetCDF weight file and converts it to
    pub fn from_weight_file(path: impl AsRef<Path> + Clone) -> Result<Self, String> {
        // open the weight file and pull its metadata out
        let weight_netcdf = open_netcdf(&path)?;
        let json_data = json_data_from_netcdf(&weight_netcdf)?;
        let grid = read_grid(&weight_netcdf)?;

        // then collect every polyid's gridpoints
        let mut polyid_gridpoints: Vec<PolyidEntry> = Vec::with_capacity(json_data.polyids.len());
        for_each_polyid(
            &weight_netcdf,
            &grid,
            json_data.polyids.len(),
            |points| {
                polyid_gridpoints.push(PolyidEntry { data: points.to_vec() });
                Ok(())
            }
        )?;

        // `from_parts` builds the lookup table out of what we just collected
        Self::from_parts(json_data, grid.lat_len, grid.lon_len, polyid_gridpoints)
    }

    /// Converts a NetCDF weight file straight to a `.nwt` file, holding no more
    /// than one polyid's worth of grid data in memory at a time. Prefer this
    /// over [`Self::from_weight_file`] followed by [`Self::serialize_to_file`]
    /// for files too large to fit in memory
    pub fn convert_weight_file(
        input: impl AsRef<Path> + Clone,
        output: impl AsRef<Path>,
    ) -> Result<(), String> {
        let weight_netcdf = open_netcdf(&input)?;
        let json_data = json_data_from_netcdf(&weight_netcdf)?;
        let grid = read_grid(&weight_netcdf)?;

        // the header and the space the lookup table needs are both known from
        // the metadata alone, so the weights can go straight to disk
        let mut writer = NwtWriter::create(output, &json_data, grid.lat_len, grid.lon_len)?;
        for_each_polyid(
            &weight_netcdf,
            &grid,
            json_data.polyids.len(),
            |points| writer.write_points(points)
        )?;

        writer.finish()
    }

    /// Creates a new structure from its constituent parts. The lookup table is
    /// derived from `polyid_gridpoints`, so the offsets it holds always match
    /// the data they describe
    pub fn from_parts(
        json_data: JsonData,
        lat_len: u64,
        lon_len: u64,
        polyid_gridpoints: Vec<PolyidEntry>
    ) -> Result<Self, String> {
        // `serialize_to_file` writes the polyid count out of the JSON data while
        // `from_nwt` uses that same count to size the lookup table it reads. If
        // the two disagree the resulting file cannot be read back at all, so
        // reject the mismatch here rather than writing a broken file
        if json_data.polyids.len() != polyid_gridpoints.len() {
            return Err(format!(
                "Mismatched parts: {} polyids but {} gridpoint entries",
                json_data.polyids.len(),
                polyid_gridpoints.len()
            ));
        }

        // build the lookup table from the data we were handed
        let mut lookup_table: Vec<(u64, u64)> = Vec::with_capacity(polyid_gridpoints.len());
        let mut running_total: u64 = 0;
        for entry in polyid_gridpoints.iter() {
            let entry_size = entry.data.len() as u64;
            lookup_table.push((running_total, entry_size));
            running_total += entry_size;
        }

        Ok(Self {
            json_data,
            lat_len,
            lon_len,
            polyid_gridpoints,
            lookup_table
        })
    }

    /// create new structure from .NWT file
    pub fn from_nwt(path: impl AsRef<Path> + Clone) -> Result<Self, String> {
        let mut reader = NwtReader::open(path)?;

        let mut polyid_gridpoints: Vec<PolyidEntry> = Vec::with_capacity(reader.polyid_count());
        for index in 0..reader.polyid_count() {
            polyid_gridpoints.push(reader.read_entry(index)?);
        }

        Ok(Self {
            lat_len: reader.lat_len,
            lon_len: reader.lon_len,
            lookup_table: reader.lookup_table.clone(),
            json_data: reader.into_json_data(),
            polyid_gridpoints
        })
    }

    /// Generically opens a weight file. If it is a NetCDF file, it is converted
    /// to the NWT format. Otherwise it is opened as standard
    pub fn open(path: impl AsRef<Path> + Clone) -> Result<Self, String> {
        if is_nwt(&path)? {
            return Self::from_nwt(path);
        }

        let new_path = PathBuf::from(format!("{}.nwt", path.as_ref().display()));
        println!(
            "[libNextWeightFile] Serializing new weight file to {}. Use this next time to avoid precomputation step",
            new_path.display()
        );
        Self::convert_weight_file(path, &new_path)?;
        Self::from_nwt(new_path)
    }

    /// Returns a dummy weight file
    pub fn dummy(input_file: impl AsRef<Path> + Clone) -> Result<Self, String> {
        let weight_netcdf = open_netcdf(&input_file)?;

        Ok(Self {
            json_data: json_data_from_netcdf(&weight_netcdf)?,
            lat_len: 0,
            lon_len: 0,
            polyid_gridpoints: Vec::new(),
            lookup_table: Vec::new()
        })
    }

    /// serializes the new weight file to disk
    pub fn serialize_to_file(&self, filename: Option<String>) -> Result<(), String> {
        // first determine our filename. Default is "test.nwt"
        let fname = match filename {
            Some(a) => a,
            None => "test.nwt".to_string()
        };

        let mut writer = NwtWriter::create(&fname, &self.json_data, self.lat_len, self.lon_len)?;
        for entry in self.polyid_gridpoints.iter() {
            writer.write_points(&entry.data)?;
        }

        writer.finish()
    }

    /// Returns the JSON metadata associated with the file
    pub fn get_json_data(&self) -> &JsonData {
        &self.json_data
    }

    /// Returns all global attributes in the file
    pub fn get_global_attrs(&self) -> &Vec<(String, String)> {
        &self.json_data.global_attrs
    }

    /// Returns all attributes associated with a given variable
    pub fn get_var_attrs(&self, var: String) -> Option<&Vec<(String,String)>> {
        self.json_data.per_variable_attrs.get(&var)
    }

    /// Returns a list of polyids
    pub fn get_polyids(&self) -> &Vec<String> {
        &self.json_data.polyids
    }

    /// Returns a reference to all grid points in the weight file
    pub fn get_gridpoints(&self) -> &Vec<PolyidEntry> {
        &self.polyid_gridpoints
    }

    /// Returns a reference to the data lookup table
    pub fn get_lookup_table(&self) -> &Vec<(u64,u64)> {
        &self.lookup_table
    }

    /// Returns the dimensions of the weight file
    pub fn get_dimensions(&self) -> (u64, u64) {
        (self.lat_len, self.lon_len)
    }

    /// Returns a raw representation of gridpoints.
    pub fn get_raw_gridpoints(&self) -> Vec<GridPoint> {
        let mut ret = Vec::new();

        // loop over each gridpoint
        for region in self.polyid_gridpoints.iter() {
            for point in region.data.iter() {
                ret.push((
                    point.0,
                    point.1,
                    point.2,
                    point.3,
                    point.4
                ))
            }
        }

        ret
    }
}

/// Writes a weight file a polyid at a time, so that a file larger than memory
/// can be produced without ever holding all of it
///
/// The metadata fixes both the header and the size of the lookup table, so the
/// gridpoints can be appended as they arrive. [`Self::finish`] goes back and
/// fills the lookup table in once every entry's size is known, and must be
/// called or the file is left without one
#[derive(Debug)]
pub struct NwtWriter {
    file: BufWriter<File>,
    path: PathBuf,
    /// where the lookup table sits, for `finish` to seek back to
    lookup_offset: u64,
    /// how many polyids the metadata promises
    expected: usize,
    lookup_table: Vec<(u64, u64)>,
    /// how many gridpoints have been written so far
    running_total: u64,
    /// reused between writes so each entry does not allocate
    buffer: Vec<u8>,
}

impl NwtWriter {
    /// Creates a weight file and writes everything that the metadata alone
    /// determines: the header, the attributes, and room for the lookup table
    pub fn create(
        path: impl AsRef<Path>,
        json_data: &JsonData,
        lat_len: u64,
        lon_len: u64,
    ) -> Result<Self, String> {
        let path = path.as_ref().to_path_buf();
        let serialized = match serde_json::to_string(json_data) {
            Ok(a) => a,
            Err(e) => return Err(format!("Failed to serialize weight file metadata: {}", e))
        };
        let num_polyids = json_data.polyids.len();

        let file = match File::create(&path) {
            Ok(a) => a,
            Err(e) => return Err(format!("Failed to open file for serialization: {}", e))
        };
        let mut file = BufWriter::with_capacity(STREAM_BUFFER_LEN, file);

        let json_offset = HEADER_LEN as u64;
        let lookup_offset = json_offset + serialized.len() as u64;

        let mut header: Vec<u8> = Vec::with_capacity(HEADER_LEN);
        header.extend_from_slice(MAGIC);
        header.extend_from_slice(&(serialized.len() as u64).to_le_bytes());
        header.extend_from_slice(&(num_polyids as u64).to_le_bytes());
        header.extend_from_slice(&lat_len.to_le_bytes());
        header.extend_from_slice(&lon_len.to_le_bytes());
        header.extend_from_slice(&json_offset.to_le_bytes());
        header.extend_from_slice(&lookup_offset.to_le_bytes());

        write_all(&mut file, &header, &path)?;
        write_all(&mut file, serialized.as_bytes(), &path)?;
        // the lookup table can only be filled in once the entries have been
        // written, so leave the room it needs
        write_all(&mut file, &vec![0u8; num_polyids * LOOKUP_ENTRY_LEN], &path)?;

        Ok(Self {
            file,
            path,
            lookup_offset,
            expected: num_polyids,
            lookup_table: Vec::with_capacity(num_polyids),
            running_total: 0,
            buffer: Vec::new()
        })
    }

    /// Appends one polyid's gridpoints
    pub fn write_entry(&mut self, entry: &PolyidEntry) -> Result<(), String> {
        self.write_points(&entry.data)
    }

    /// Appends one polyid's gridpoints
    pub fn write_points(&mut self, points: &[GridPoint]) -> Result<(), String> {
        if self.lookup_table.len() >= self.expected {
            return Err(format!(
                "Tried to write more than the {} polyids the metadata declares",
                self.expected
            ));
        }

        self.buffer.clear();
        self.buffer.reserve(points.len() * RECORD_LEN);
        for point in points.iter() {
            self.buffer.extend_from_slice(&point.0.to_le_bytes());
            self.buffer.extend_from_slice(&point.1.to_le_bytes());
            self.buffer.extend_from_slice(&point.2.to_le_bytes());
            self.buffer.extend_from_slice(&point.3.to_le_bytes());
            self.buffer.extend_from_slice(&point.4.to_le_bytes());
        }
        write_all(&mut self.file, &self.buffer, &self.path)?;

        self.lookup_table.push((self.running_total, points.len() as u64));
        self.running_total += points.len() as u64;

        Ok(())
    }

    /// How many polyids have been written so far
    pub fn written(&self) -> usize {
        self.lookup_table.len()
    }

    /// Fills in the lookup table and closes the file
    pub fn finish(mut self) -> Result<(), String> {
        if self.lookup_table.len() != self.expected {
            return Err(format!(
                "Wrote {} polyids but the metadata declares {}",
                self.lookup_table.len(),
                self.expected
            ));
        }

        let mut table: Vec<u8> = Vec::with_capacity(self.lookup_table.len() * LOOKUP_ENTRY_LEN);
        for (offset, count) in self.lookup_table.iter() {
            table.extend_from_slice(&offset.to_le_bytes());
            table.extend_from_slice(&count.to_le_bytes());
        }

        // go back to the room left for the table in `create`
        if let Err(e) = self.file.seek(SeekFrom::Start(self.lookup_offset)) {
            return Err(format!("Failed to seek in {}: {}", self.path.display(), e));
        }
        write_all(&mut self.file, &table, &self.path)?;

        if let Err(e) = self.file.flush() {
            return Err(format!("Failed to flush {}: {}", self.path.display(), e));
        }

        Ok(())
    }
}

/// Reads a weight file a polyid at a time, so that a file larger than memory
/// can be worked through without ever holding all of it
///
/// The header, the attributes and the lookup table are read up front, since
/// they are small and the table is what makes the gridpoints addressable.
/// Everything else is read on demand
#[derive(Debug)]
pub struct NwtReader {
    file: BufReader<File>,
    path: PathBuf,
    json_data: JsonData,
    lat_len: u64,
    lon_len: u64,
    lookup_table: Vec<(u64, u64)>,
    /// where the gridpoints start, just past the lookup table
    data_offset: u64,
    /// where the file is positioned, so sequential reads can skip the seek and
    /// let the buffer do its job
    position: u64,
    /// reused between reads so each entry does not allocate
    buffer: Vec<u8>,
}

impl NwtReader {
    /// Opens a weight file, reading its metadata but none of its gridpoints
    pub fn open(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref().to_path_buf();
        let file = match File::open(&path) {
            Ok(a) => a,
            Err(e) => return Err(format!("Failed to open {}: {}", path.display(), e))
        };
        let file_len = match file.metadata() {
            Ok(a) => a.len(),
            Err(e) => return Err(format!("Failed to stat {}: {}", path.display(), e))
        };
        let mut file = BufReader::with_capacity(STREAM_BUFFER_LEN, file);

        // the magic bytes first: a file that does not start with them is not a
        // weight file at all, which is worth saying before anything about its
        // length or contents
        let mut magic = [0u8; 4];
        match file.read_exact(&mut magic) {
            Ok(()) if &magic == MAGIC => {}
            Ok(()) => return Err("Invalid file format".to_string()),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                return Err("Invalid file format".to_string())
            }
            Err(e) => return Err(format!("Failed to read {}: {}", path.display(), e))
        }

        // then the rest of the header, which says where everything else lives
        let mut header = [0u8; HEADER_LEN - 4];
        read_exact(&mut file, &mut header, &path, "header")?;

        let json_len = le_u64(&header, 0);
        let num_polyids = le_u64(&header, 8);
        let lat_len = le_u64(&header, 16);
        let lon_len = le_u64(&header, 24);
        let json_offset = le_u64(&header, 32);
        let lookup_offset = le_u64(&header, 40);

        // a corrupt header can claim any size at all, so check the offsets
        // against the file before trusting them with an allocation
        let lookup_len = num_polyids * LOOKUP_ENTRY_LEN as u64;
        check_extent(json_offset, json_len, file_len, &path, "attributes")?;
        check_extent(lookup_offset, lookup_len, file_len, &path, "lookup table")?;

        // then the attributes...
        let mut json_bytes = vec![0u8; json_len as usize];
        seek(&mut file, json_offset, &path)?;
        read_exact(&mut file, &mut json_bytes, &path, "attributes")?;
        let json_text = match String::from_utf8(json_bytes) {
            Ok(a) => a,
            Err(e) => return Err(format!("{} holds invalid attribute text: {}", path.display(), e))
        };
        let json_data: JsonData = match serde_json::from_str(&json_text) {
            Ok(a) => a,
            Err(e) => return Err(format!("{} holds invalid attributes: {}", path.display(), e))
        };

        // ... and the lookup table, which makes the gridpoints addressable
        let mut lookup_bytes = vec![0u8; lookup_len as usize];
        seek(&mut file, lookup_offset, &path)?;
        read_exact(&mut file, &mut lookup_bytes, &path, "lookup table")?;
        let mut lookup_table: Vec<(u64, u64)> = Vec::with_capacity(num_polyids as usize);
        for entry in lookup_bytes.as_chunks::<LOOKUP_ENTRY_LEN>().0 {
            lookup_table.push((le_u64(entry, 0), le_u64(entry, 8)));
        }

        let data_offset = lookup_offset + lookup_len;

        Ok(Self {
            file,
            path,
            json_data,
            lat_len,
            lon_len,
            lookup_table,
            data_offset,
            position: data_offset,
            buffer: Vec::new()
        })
    }

    /// Reads the gridpoints of the polyid at `index`
    pub fn read_entry(&mut self, index: usize) -> Result<PolyidEntry, String> {
        let (offset, count) = match self.lookup_table.get(index) {
            Some(a) => *a,
            None => return Err(format!(
                "{} holds {} polyids, so there is none at index {}",
                self.path.display(),
                self.lookup_table.len(),
                index
            ))
        };

        let start = self.data_offset + offset * RECORD_LEN as u64;
        let len = count as usize * RECORD_LEN;

        // walking the file in order is the common case, and seeking would throw
        // the read buffer away every time
        if start != self.position {
            seek(&mut self.file, start, &self.path)?;
            self.position = start;
        }

        self.buffer.resize(len, 0);
        read_exact(&mut self.file, &mut self.buffer, &self.path, "gridpoints")?;
        self.position += len as u64;

        let mut entry = PolyidEntry::new();
        entry.data.reserve(count as usize);
        for record in self.buffer.as_chunks::<RECORD_LEN>().0 {
            entry.add_point(
                le_u32(record, 0),
                le_u32(record, 4),
                le_f32(record, 8),
                le_f32(record, 12),
                le_f32(record, 16)
            );
        }

        Ok(entry)
    }

    /// Walks every polyid in the file in order
    pub fn entries(&mut self) -> Entries<'_> {
        Entries { reader: self, index: 0 }
    }

    /// How many polyids the file holds
    pub fn polyid_count(&self) -> usize {
        self.lookup_table.len()
    }

    /// Whether the file holds no polyids at all
    pub fn is_empty(&self) -> bool {
        self.lookup_table.is_empty()
    }

    /// Returns the JSON metadata associated with the file
    pub fn get_json_data(&self) -> &JsonData {
        &self.json_data
    }

    /// Returns all global attributes in the file
    pub fn get_global_attrs(&self) -> &Vec<(String, String)> {
        &self.json_data.global_attrs
    }

    /// Returns a list of polyids
    pub fn get_polyids(&self) -> &Vec<String> {
        &self.json_data.polyids
    }

    /// Returns a reference to the data lookup table
    pub fn get_lookup_table(&self) -> &Vec<(u64, u64)> {
        &self.lookup_table
    }

    /// Returns the dimensions of the weight file
    pub fn get_dimensions(&self) -> (u64, u64) {
        (self.lat_len, self.lon_len)
    }

    /// Takes the metadata, leaving the file behind
    pub fn into_json_data(self) -> JsonData {
        self.json_data
    }
}

/// Walks the polyids of a weight file in order. Created by
/// [`NwtReader::entries`]
#[derive(Debug)]
pub struct Entries<'a> {
    reader: &'a mut NwtReader,
    index: usize,
}

impl Iterator for Entries<'_> {
    type Item = Result<PolyidEntry, String>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.reader.polyid_count() {
            return None;
        }

        let entry = self.reader.read_entry(self.index);
        self.index += 1;

        Some(entry)
    }
}

impl JsonData {
    /// creates a new instance of `JsonData`
    pub fn new() -> Self {
        Self {
            global_attrs: Vec::new(),
            per_variable_attrs: BTreeMap::new(),
            polyids: Vec::new()
        }
    }

    /// adds a global attribute to the structure
    pub fn add_global_attr(&mut self, key: String, value: String) {
        self.global_attrs.push((key, value));
    }

    /// adds a new variable to the structure
    pub fn add_variable(&mut self, variable_name: &str) {
        self.per_variable_attrs.insert(variable_name.to_string(), Vec::new());
    }

    /// adds a new attribute for the associated variable. If the variable has
    /// not yet been added, it is added
    pub fn add_variable_attr(&mut self, var_name: &str, key: String, value: String) {
        // check if it exists, if not add it
        if !self.per_variable_attrs.contains_key(var_name) {
            self.add_variable(var_name);
        }
        // add the values to the vector
        let vec_ref = self.per_variable_attrs.get_mut(var_name).unwrap();
        vec_ref.push((key, value));
    }

    /// adds a polyid to the list of polyids
    pub fn add_polyid(&mut self, polyid: String) {
        self.polyids.push(polyid)
    }

    /// Retrieves a global arribute with the provided name
    pub fn get_global_attr(&self, name: &str) -> Result<String, String> {
        for v in self.global_attrs.iter() {
            if v.0 == *name { return Ok(v.1.clone()) }
        }
        Err(format!("Global attribute {} not found", name))
    }

    /// Retrieves all global attributes
    pub fn get_global_attrs(&self) -> &Vec<(String, String)> {
        &self.global_attrs
    }

    /// Retrieves the names of every variable in the structure
    pub fn get_variables(&self) -> Vec<&String> {
        self.per_variable_attrs.keys().collect()
    }

    /// Retrieves all attributes associated with a given variable
    pub fn get_var_attrs(&self, variable_name: &str) -> Option<&Vec<(String, String)>> {
        self.per_variable_attrs.get(variable_name)
    }

    /// Retrieves a given variable's attribute of a provided name
    pub fn get_var_attr(&self, variable_name: &str, attr_name: &str) -> Result<String, String> {
        match self.per_variable_attrs.get(variable_name) {
            Some(a) => {
                for v in a.iter() {
                    if v.0 == *attr_name { return Ok(v.1.clone()) }
                }
                Err(format!("Global attribute {} not found", attr_name))
            },
            None => Err(format!("Variable {} not found in the weight file", variable_name))
        }
    }
}

impl Default for JsonData {
    fn default() -> Self {
        Self::new()
    }
}

impl PolyidEntry {
    /// creates a new PolyidEntry
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    /// adds a new point to the PolyidEntry
    pub fn add_point(&mut self, lat_idx: u32, lon_idx: u32, lat: f32, lon: f32, value: f32) {
        self.data.push((lat_idx, lon_idx, lat, lon, value))
    }
}

impl Default for PolyidEntry {
    fn default() -> Self {
        Self::new()
    }
}

/// The lat/lon grid a weight file's indices are relative to
struct Grid {
    lat_len: u64,
    lon_len: u64,
    lat_vals: Vec<f32>,
    lon_vals: Vec<f32>,
}

/// Opens a NetCDF weight file
fn open_netcdf(path: &(impl AsRef<Path> + Clone)) -> Result<netcdf::File, String> {
    match netcdf::open(path.clone()) {
        Ok(a) => Ok(a),
        Err(e) => Err(format!("Failed to open {}: {}", path.as_ref().display(), e))
    }
}

/// Fetches a variable, naming it if it is missing
fn netcdf_variable<'a>(
    weight_netcdf: &'a netcdf::File,
    name: &str,
) -> Result<netcdf::Variable<'a>, String> {
    match weight_netcdf.variable(name) {
        Some(a) => Ok(a),
        None => Err(format!("Weight file has no {} variable", name))
    }
}

/// Flattens a NetCDF attribute into the string the weight file stores
fn attribute_to_string(value: AttributeValue) -> String {
    match value {
        AttributeValue::Str(a) => a,
        AttributeValue::Strs(a) => a[0].clone(),
        AttributeValue::Double(a) => format!("{}", a),
        AttributeValue::Doubles(a) => format!("{:?}", a),
        AttributeValue::Float(a) => format!("{}", a),
        AttributeValue::Floats(a) => format!("{:?}", a),
        AttributeValue::Int(a) => format!("{}", a),
        AttributeValue::Ints(a) => format!("{:?}", a),
        AttributeValue::Longlong(a) => format!("{}", a),
        AttributeValue::Longlongs(a) => format!("{:?}", a),
        AttributeValue::Schar(a) => format!("{}", a),
        AttributeValue::Schars(a) => format!("{:?}", a),
        AttributeValue::Short(a) => format!("{}", a),
        AttributeValue::Shorts(a) => format!("{:?}", a),
        AttributeValue::Uchar(a) => format!("{}", a),
        AttributeValue::Uchars(a) => format!("{:?}", a),
        AttributeValue::Uint(a) => format!("{}", a),
        AttributeValue::Uints(a) => format!("{:?}", a),
        AttributeValue::Ulonglong(a) => format!("{}", a),
        AttributeValue::Ulonglongs(a) => format!("{:?}", a),
        AttributeValue::Ushort(a) => format!("{}", a),
        AttributeValue::Ushorts(a) => format!("{:?}", a)
    }
}

/// Collects the attributes and polyid names of a NetCDF weight file
fn json_data_from_netcdf(weight_netcdf: &netcdf::File) -> Result<JsonData, String> {
    let mut json_data = JsonData::new();

    // now we get all of the attributes...
    for attr in weight_netcdf.attributes() {
        let attr_value = match attr.value() {
            Ok(a) => attribute_to_string(a),
            Err(e) => return Err(format!("Failed to read global attribute {}: {}", attr.name(), e))
        };
        // add it to our list of global attributes
        json_data.add_global_attr(attr.name().to_string(), attr_value);
    }

    // ... and add all of the variable attributes ...
    for var in weight_netcdf.variables() {
        let var_name = var.name();
        json_data.add_variable(&var_name);
        for attr in var.attributes() {
            if attr.name() == "_FillValue" {
                continue;
            }
            let attr_value = match attr.value() {
                Ok(a) => attribute_to_string(a),
                Err(e) => return Err(format!(
                    "Failed to read attribute {} of variable {}: {}",
                    attr.name(), var_name, e
                ))
            };
            json_data.add_variable_attr(&var_name, attr.name().to_string(), attr_value);
        }
    }

    // now that we have gotten our attributes all squared away, lets start
    // looking at data. First things first, lets store those polyids
    let polyid_var = netcdf_variable(weight_netcdf, "polyid")?;
    for polyid in 0..polyid_var.len() {
        match polyid_var.get_string(polyid) {
            Ok(a) => json_data.add_polyid(a),
            Err(e) => return Err(format!("Failed to read polyid {}: {}", polyid, e))
        }
    }

    Ok(json_data)
}

/// Reads the lat/lon grid of a NetCDF weight file
fn read_grid(weight_netcdf: &netcdf::File) -> Result<Grid, String> {
    let latvar = netcdf_variable(weight_netcdf, "lat")?;
    let lonvar = netcdf_variable(weight_netcdf, "lon")?;

    let lat_vals = match latvar.get_values::<f32, _>(..) {
        Ok(a) => a,
        Err(e) => return Err(format!("Failed to read latitudes: {}", e))
    };
    let lon_vals = match lonvar.get_values::<f32, _>(..) {
        Ok(a) => a,
        Err(e) => return Err(format!("Failed to read longitudes: {}", e))
    };

    let lat_len = match weight_netcdf.dimension("lat") {
        Some(a) => a.len() as u64,
        None => return Err("Weight file has no lat dimension".to_string())
    };
    let lon_len = match weight_netcdf.dimension("lon") {
        Some(a) => a.len() as u64,
        None => return Err("Weight file has no lon dimension".to_string())
    };

    Ok(Grid { lat_len, lon_len, lat_vals, lon_vals })
}

/// Hands every polyid's gridpoints to `handle`, one polyid at a time. Only the
/// weights of the polyid in hand are held in memory, which is what lets a file
/// bigger than memory be converted
fn for_each_polyid<F>(
    weight_netcdf: &netcdf::File,
    grid: &Grid,
    num_polyids: usize,
    mut handle: F,
) -> Result<(), String>
where
    F: FnMut(&[GridPoint]) -> Result<(), String>,
{
    let regridweights = netcdf_variable(weight_netcdf, "regridweights")?;
    // a file without a fill value has nothing to leave out
    let fill = match regridweights.fill_value::<f32>() {
        Ok(a) => a,
        Err(e) => return Err(format!("Failed to read the regridweights fill value: {}", e))
    };

    let lat_len = grid.lat_len as usize;
    let lon_len = grid.lon_len as usize;
    let mut points: Vec<GridPoint> = Vec::new();

    // for every polyid...
    for polyid in 0..num_polyids {
        let data = match regridweights.get::<f32, _>((polyid, .., ..)) {
            Ok(a) => a,
            Err(e) => return Err(format!("Failed to read the weights of polyid {}: {}", polyid, e))
        };
        let dat_slice = match data.as_slice() {
            Some(a) => a,
            None => return Err(format!("The weights of polyid {} are not contiguous", polyid))
        };

        points.clear();
        // ... for every data value...
        for lat_idx in 0..lat_len {
            for lon_idx in 0..lon_len {
                let data_value = dat_slice[lat_idx * lon_len + lon_idx];
                // ...if it isnt a fill value...
                if Some(data_value) != fill {
                    // ... then calculate the lat lon and save the weight
                    points.push((
                        lat_idx as u32,
                        lon_idx as u32,
                        grid.lat_vals[lat_idx],
                        grid.lon_vals[lon_idx],
                        data_value
                    ));
                }
            }
        }

        handle(&points)?;
    }

    Ok(())
}

/// Checks a file for the weight file magic bytes
fn is_nwt(path: &(impl AsRef<Path> + Clone)) -> Result<bool, String> {
    let path = path.as_ref();
    let mut magic = [0u8; 4];

    // scope brackets here to make sure the file is closed before it is reopened
    let mut input_file = match File::open(path) {
        Ok(a) => a,
        Err(e) => return Err(format!("Failed to open {}: {}", path.display(), e))
    };
    match input_file.read_exact(&mut magic) {
        Ok(()) => Ok(&magic == MAGIC),
        // too small to hold a header, so it is not one of ours
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(format!("Failed to read {}: {}", path.display(), e))
    }
}

/// Writes a whole buffer, naming the file if it fails
fn write_all(file: &mut BufWriter<File>, bytes: &[u8], path: &Path) -> Result<(), String> {
    match file.write_all(bytes) {
        Ok(()) => Ok(()),
        Err(e) => Err(format!("Failed to write to {}: {}", path.display(), e))
    }
}

/// Fills a whole buffer, naming the file if it fails
fn read_exact(
    file: &mut BufReader<File>,
    bytes: &mut [u8],
    path: &Path,
    what: &str,
) -> Result<(), String> {
    match file.read_exact(bytes) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Err(format!(
            "{} is corrupt or truncated: it ends inside the {}",
            path.display(), what
        )),
        Err(e) => Err(format!("Failed to read the {} of {}: {}", what, path.display(), e))
    }
}

/// Moves to a position in the file, naming it if it fails
fn seek(file: &mut BufReader<File>, to: u64, path: &Path) -> Result<(), String> {
    match file.seek(SeekFrom::Start(to)) {
        Ok(_) => Ok(()),
        Err(e) => Err(format!("Failed to seek in {}: {}", path.display(), e))
    }
}

/// Checks that a stretch of the file the header points at is really there. A
/// corrupt header can claim any size at all, and sizing a buffer from one
/// would ask for an absurd allocation
fn check_extent(
    offset: u64,
    len: u64,
    file_len: u64,
    path: &Path,
    what: &str,
) -> Result<(), String> {
    if offset.saturating_add(len) > file_len {
        return Err(format!(
            "{} is corrupt or truncated: its {} needs {} bytes at offset {} but the file is {} bytes",
            path.display(), what, len, offset, file_len
        ));
    }

    Ok(())
}

/// Reads a little endian `u64` out of a buffer
fn le_u64(bytes: &[u8], at: usize) -> u64 {
    let mut buff = [0u8; size_of::<u64>()];
    buff.copy_from_slice(&bytes[at..at + size_of::<u64>()]);
    u64::from_le_bytes(buff)
}

/// Reads a little endian `u32` out of a buffer
fn le_u32(bytes: &[u8], at: usize) -> u32 {
    let mut buff = [0u8; size_of::<u32>()];
    buff.copy_from_slice(&bytes[at..at + size_of::<u32>()]);
    u32::from_le_bytes(buff)
}

/// Reads a little endian `f32` out of a buffer
fn le_f32(bytes: &[u8], at: usize) -> f32 {
    let mut buff = [0u8; size_of::<f32>()];
    buff.copy_from_slice(&bytes[at..at + size_of::<f32>()]);
    f32::from_le_bytes(buff)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The fill value the NetCDF fixture leaves in cells with no weight
    const FILL: f32 = -9999.0;

    /// A scratch path that will not collide with another test's
    fn scratch(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!("nextweight_test_{}", name));
        let _ = std::fs::remove_file(&path);
        path
    }

    /// A small weight file built in memory
    fn sample(polyids: &[(&str, Vec<GridPoint>)]) -> NextWeightFile {
        let mut json_data = JsonData::new();
        json_data.add_global_attr("title".to_string(), "sample".to_string());
        json_data.add_variable_attr("regridweights", "units".to_string(), "1".to_string());

        let mut entries = Vec::new();
        for (polyid, points) in polyids.iter() {
            json_data.add_polyid(polyid.to_string());
            entries.push(PolyidEntry { data: points.clone() });
        }

        NextWeightFile::from_parts(json_data, 180, 360, entries).unwrap()
    }

    #[test]
    fn json_header_is_independent_of_insertion_order() {
        let vars = ["regridweights", "lat", "lon", "lat_bnds", "lon_bnds", "polyid", "bnds"];

        let build = |order: &mut dyn Iterator<Item = &&str>| {
            let mut json_data = JsonData::new();
            for var in order {
                json_data.add_variable_attr(var, "long_name".to_string(), var.to_string());
            }
            serde_json::to_string(&json_data).unwrap()
        };

        assert_eq!(build(&mut vars.iter()), build(&mut vars.iter().rev()));
    }

    #[test]
    #[ignore = "needs the weight files under ../rust_science/test_cases"]
    fn it_works() {
        // lets test this
        let test_path = PathBuf::from("../rust_science/test_cases/wgts/cckp_aggregation_1x1");
        let new_path = PathBuf::from("test.nwt");
        let new_weight = NextWeightFile::from_weight_file(test_path).unwrap();
        println!("new_weight file weight examples: {:?}", &new_weight.lookup_table[..10]);
        println!("new_weight file weights: {:?}", &new_weight.polyid_gridpoints[..10]);
        new_weight.serialize_to_file(None).unwrap();
        let fresh_weight = NextWeightFile::from_nwt(new_path).unwrap();

        for v in 0..new_weight.json_data.polyids.len() {
            assert_eq!(new_weight.json_data.polyids[v], fresh_weight.json_data.polyids[v]);
        }

        for v in 0..new_weight.json_data.global_attrs.len() {
            let attr_name = &new_weight.json_data.global_attrs[v].0;
            let attr_value = &new_weight.json_data.global_attrs[v].1;
            let new_name = &fresh_weight.json_data.global_attrs[v].0;
            let new_value = &fresh_weight.json_data.global_attrs[v].1;

            assert_eq!(attr_name, new_name);
            assert_eq!(attr_value, new_value);
        }

        for key in new_weight.json_data.per_variable_attrs.keys() {
            let curr_v_vec = new_weight.json_data.per_variable_attrs.get(key).unwrap();
            let new_v_vec = fresh_weight.json_data.per_variable_attrs.get(key).unwrap();
            for v in 0..curr_v_vec.len() {
                let attr_name = &curr_v_vec[v].0;
                let attr_value = &curr_v_vec[v].1;
                let new_name = &new_v_vec[v].0;
                let new_value = &new_v_vec[v].1;

                assert_eq!(attr_name, new_name);
                assert_eq!(attr_value, new_value);
            }
        }

        for v in 0..new_weight.lookup_table.len() {
            let curr_v = new_weight.lookup_table[v];
            let new_v = fresh_weight.lookup_table[v];
            assert_eq!(curr_v, new_v);
        }

        for v in 0..new_weight.polyid_gridpoints.len() {
            let curr_v_entry = &new_weight.polyid_gridpoints[v];
            let new_v_entry = &fresh_weight.polyid_gridpoints[v];
            for grid_idx in 0..curr_v_entry.data.len() {
                let curr_pt = curr_v_entry.data[grid_idx];
                let new_pt = new_v_entry.data[grid_idx];

                let curr_lat = curr_pt.0;
                let curr_lon = curr_pt.1;
                let curr_wgt = curr_pt.2;
                let new_lon = new_pt.1;
                let new_lat = new_pt.0;
                let new_wgt = new_pt.2;

                assert_eq!(curr_lat, new_lat);
                assert_eq!(curr_lon, new_lon);
                assert_eq!(curr_wgt, new_wgt);
            }
        }

    }

    #[test]
    fn from_parts_round_trip() {
        let out_path = scratch("from_parts.nwt");

        // build a small weight file by hand...
        let built = sample(&[
            ("AAA", vec![(0, 1, -89.5, -179.5, 0.25), (1, 2, -88.5, -178.5, 0.75)]),
            ("BBB", vec![(3, 4, -86.5, -176.5, 1.0)])
        ]);

        // ... the derived lookup table should describe the data that follows it ...
        assert_eq!(built.get_lookup_table(), &vec![(0u64, 2u64), (2u64, 1u64)]);
        assert_eq!(built.get_json_data().get_variables(), vec!["regridweights"]);

        // ... and it should survive a trip through the on-disk format
        built.serialize_to_file(Some(out_path.to_str().unwrap().to_string())).unwrap();
        let read_back = NextWeightFile::from_nwt(&out_path).unwrap();

        assert_eq!(read_back.get_dimensions(), (180, 360));
        assert_eq!(read_back.get_polyids(), built.get_polyids());
        assert_eq!(read_back.get_lookup_table(), built.get_lookup_table());
        assert_eq!(read_back.get_raw_gridpoints(), built.get_raw_gridpoints());
        assert_eq!(
            read_back.get_json_data().get_var_attrs("regridweights"),
            built.get_json_data().get_var_attrs("regridweights")
        );

        std::fs::remove_file(out_path).unwrap();
    }

    #[test]
    fn from_parts_rejects_mismatched_parts() {
        let mut json_data = JsonData::new();
        json_data.add_polyid("AAA".to_string());

        // a polyid with no gridpoint entry to go with it would serialize to a
        // file that cannot be read back
        assert!(NextWeightFile::from_parts(json_data, 1, 1, Vec::new()).is_err());
    }

    /// Pins the on-disk layout. Anything that shifts a field breaks every file
    /// already written, so the bytes are checked against the format by hand
    #[test]
    fn writes_the_documented_layout() {
        let out_path = scratch("layout.nwt");

        let mut json_data = JsonData::new();
        json_data.add_polyid("AAA".to_string());
        let serialized = serde_json::to_string(&json_data).unwrap();

        let mut writer = NwtWriter::create(&out_path, &json_data, 7, 9).unwrap();
        writer.write_points(&[(1, 2, 3.5, 4.5, 5.5)]).unwrap();
        assert_eq!(writer.written(), 1);
        writer.finish().unwrap();

        let bytes = std::fs::read(&out_path).unwrap();
        let json_offset = 4 + 6 * 8;
        let lookup_offset = json_offset + serialized.len();

        assert_eq!(&bytes[0..4], b"NEWT");
        assert_eq!(le_u64(&bytes, 4), serialized.len() as u64, "json length");
        assert_eq!(le_u64(&bytes, 12), 1, "polyid count");
        assert_eq!(le_u64(&bytes, 20), 7, "lat length");
        assert_eq!(le_u64(&bytes, 28), 9, "lon length");
        assert_eq!(le_u64(&bytes, 36), json_offset as u64, "json offset");
        assert_eq!(le_u64(&bytes, 44), lookup_offset as u64, "lookup offset");
        assert_eq!(
            &bytes[json_offset..lookup_offset],
            serialized.as_bytes(),
            "attributes"
        );

        // one lookup entry, then the point it points at
        assert_eq!(le_u64(&bytes, lookup_offset), 0, "first data offset");
        assert_eq!(le_u64(&bytes, lookup_offset + 8), 1, "first point count");
        let data = lookup_offset + 16;
        assert_eq!(bytes.len(), data + 20, "one 20 byte record");
        assert_eq!(le_u32(&bytes, data), 1);
        assert_eq!(le_u32(&bytes, data + 4), 2);
        assert_eq!(le_f32(&bytes, data + 8), 3.5);
        assert_eq!(le_f32(&bytes, data + 12), 4.5);
        assert_eq!(le_f32(&bytes, data + 16), 5.5);

        std::fs::remove_file(out_path).unwrap();
    }

    #[test]
    fn streaming_reader_matches_reading_the_whole_file() {
        let out_path = scratch("streaming.nwt");

        let built = sample(&[
            ("AAA", vec![(0, 1, -89.5, -179.5, 0.25), (1, 2, -88.5, -178.5, 0.75)]),
            ("BBB", vec![]),
            ("CCC", vec![(3, 4, -86.5, -176.5, 1.0)])
        ]);
        built.serialize_to_file(Some(out_path.to_str().unwrap().to_string())).unwrap();

        let whole = NextWeightFile::from_nwt(&out_path).unwrap();
        let mut reader = NwtReader::open(&out_path).unwrap();

        assert_eq!(reader.polyid_count(), 3);
        assert!(!reader.is_empty());
        assert_eq!(reader.get_polyids(), whole.get_polyids());
        assert_eq!(reader.get_dimensions(), whole.get_dimensions());
        assert_eq!(reader.get_lookup_table(), whole.get_lookup_table());

        // walking the file must give back exactly what loading it does, empty
        // polyids and all
        let streamed: Vec<PolyidEntry> = reader
            .entries()
            .collect::<Result<Vec<_>, String>>()
            .unwrap();
        let loaded = whole.get_gridpoints();
        assert_eq!(streamed.len(), loaded.len());
        for (streamed, loaded) in streamed.iter().zip(loaded.iter()) {
            assert_eq!(streamed.data, loaded.data);
        }

        // and so must reading one straight out of the middle
        assert_eq!(reader.read_entry(2).unwrap().data, loaded[2].data);
        assert!(reader.read_entry(3).is_err(), "there is no fourth polyid");

        std::fs::remove_file(out_path).unwrap();
    }

    #[test]
    fn reader_rejects_a_truncated_file() {
        let out_path = scratch("truncated.nwt");

        let built = sample(&[("AAA", vec![(0, 1, -89.5, -179.5, 0.25)])]);
        built.serialize_to_file(Some(out_path.to_str().unwrap().to_string())).unwrap();

        let bytes = std::fs::read(&out_path).unwrap();

        // cut off inside the header
        std::fs::write(&out_path, &bytes[..12]).unwrap();
        let err = NwtReader::open(&out_path).unwrap_err();
        assert!(err.contains("corrupt or truncated"), "{}", err);

        // cut off inside the attributes, which the header still promises
        std::fs::write(&out_path, &bytes[..bytes.len() - 25]).unwrap();
        let err = NwtReader::open(&out_path).unwrap_err();
        assert!(err.contains("corrupt or truncated"), "{}", err);

        // and something else entirely
        std::fs::write(&out_path, b"PK\x03\x04 this is a zip file").unwrap();
        assert_eq!(NwtReader::open(&out_path).unwrap_err(), "Invalid file format");

        std::fs::remove_file(out_path).unwrap();
    }

    #[test]
    fn writer_rejects_the_wrong_number_of_polyids() {
        let out_path = scratch("miscounted.nwt");

        let mut json_data = JsonData::new();
        json_data.add_polyid("AAA".to_string());
        json_data.add_polyid("BBB".to_string());

        // one polyid short of what the metadata promises
        let mut writer = NwtWriter::create(&out_path, &json_data, 1, 1).unwrap();
        writer.write_points(&[(0, 0, 0.0, 0.0, 1.0)]).unwrap();
        assert!(writer.finish().is_err());

        // and one too many
        let mut writer = NwtWriter::create(&out_path, &json_data, 1, 1).unwrap();
        writer.write_points(&[(0, 0, 0.0, 0.0, 1.0)]).unwrap();
        writer.write_points(&[(0, 0, 0.0, 0.0, 1.0)]).unwrap();
        assert!(writer.write_points(&[(0, 0, 0.0, 0.0, 1.0)]).is_err());

        std::fs::remove_file(out_path).unwrap();
    }

    /// Writes a small NetCDF weight file: 3 polyids over a 4x5 grid, where only
    /// every third cell carries a weight
    fn write_netcdf_fixture(path: &Path) -> Vec<f32> {
        let mut file = netcdf::create(path).unwrap();
        file.add_dimension("polyid", 3).unwrap();
        file.add_dimension("lat", 4).unwrap();
        file.add_dimension("lon", 5).unwrap();
        file.add_attribute("title", "netcdf fixture").unwrap();

        let mut lat = file.add_variable::<f32>("lat", &["lat"]).unwrap();
        lat.put_values(&[-67.5f32, -22.5, 22.5, 67.5], ..).unwrap();
        lat.put_attribute("units", "degrees_north").unwrap();

        let mut lon = file.add_variable::<f32>("lon", &["lon"]).unwrap();
        lon.put_values(&[-144.0f32, -72.0, 0.0, 72.0, 144.0], ..).unwrap();
        lon.put_attribute("units", "degrees_east").unwrap();

        let mut polyid = file.add_string_variable("polyid", &["polyid"]).unwrap();
        for (index, name) in ["AAA", "BBB", "CCC"].iter().enumerate() {
            polyid.put_string(name, index).unwrap();
        }

        let mut values: Vec<f32> = Vec::with_capacity(3 * 4 * 5);
        for i in 0..(3 * 4 * 5) {
            values.push(if i % 3 == 0 { 0.1 * i as f32 } else { FILL });
        }

        let mut weights = file
            .add_variable::<f32>("regridweights", &["polyid", "lat", "lon"])
            .unwrap();
        weights.set_fill_value(FILL).unwrap();
        weights.put_attribute("units", "1").unwrap();
        weights.put_values(&values, ..).unwrap();

        values
    }

    #[test]
    fn converts_a_netcdf_weight_file_without_loading_it() {
        let nc_path = scratch("fixture.nc");
        let nwt_path = scratch("converted.nwt");
        let values = write_netcdf_fixture(&nc_path);

        // the streaming conversion and the in-memory one have to agree
        NextWeightFile::convert_weight_file(&nc_path, &nwt_path).unwrap();
        let streamed = NextWeightFile::from_nwt(&nwt_path).unwrap();
        let in_memory = NextWeightFile::from_weight_file(&nc_path).unwrap();

        assert_eq!(streamed.get_polyids(), &vec![
            "AAA".to_string(),
            "BBB".to_string(),
            "CCC".to_string()
        ]);
        assert_eq!(streamed.get_dimensions(), (4, 5));
        assert_eq!(streamed.get_polyids(), in_memory.get_polyids());
        assert_eq!(streamed.get_dimensions(), in_memory.get_dimensions());
        assert_eq!(streamed.get_lookup_table(), in_memory.get_lookup_table());
        assert_eq!(streamed.get_raw_gridpoints(), in_memory.get_raw_gridpoints());
        assert_eq!(streamed.get_global_attrs(), in_memory.get_global_attrs());

        // the fill values are left out, and the rest are placed on the grid
        let points = streamed.get_raw_gridpoints();
        assert_eq!(points.len(), values.iter().filter(|v| **v != FILL).count());
        let lats = [-67.5f32, -22.5, 22.5, 67.5];
        let lons = [-144.0f32, -72.0, 0.0, 72.0, 144.0];
        for (polyid, entry) in streamed.get_gridpoints().iter().enumerate() {
            for (lat_idx, lon_idx, lat, lon, weight) in entry.data.iter() {
                assert_eq!(*lat, lats[*lat_idx as usize]);
                assert_eq!(*lon, lons[*lon_idx as usize]);
                assert_eq!(*weight, values[polyid * 20 + *lat_idx as usize * 5 + *lon_idx as usize]);
            }
        }

        // the attributes come across, bar the fill value
        assert_eq!(
            streamed.get_json_data().get_global_attr("title"),
            Ok("netcdf fixture".to_string())
        );
        assert_eq!(
            streamed.get_json_data().get_var_attrs("regridweights"),
            Some(&vec![("units".to_string(), "1".to_string())])
        );
        assert_eq!(
            streamed.get_json_data().get_var_attr("lat", "units"),
            Ok("degrees_north".to_string())
        );

        std::fs::remove_file(nc_path).unwrap();
        std::fs::remove_file(nwt_path).unwrap();
    }

    #[test]
    fn open_converts_a_netcdf_file_once() {
        let nc_path = scratch("open_fixture.nc");
        write_netcdf_fixture(&nc_path);
        let cached = PathBuf::from(format!("{}.nwt", nc_path.display()));
        let _ = std::fs::remove_file(&cached);

        // the first open converts and leaves the converted file behind...
        let first = NextWeightFile::open(&nc_path).unwrap();
        assert!(cached.exists(), "the converted file should be kept");

        // ... which the second open reads instead of converting again
        let second = NextWeightFile::open(&cached).unwrap();
        assert_eq!(first.get_polyids(), second.get_polyids());
        assert_eq!(first.get_raw_gridpoints(), second.get_raw_gridpoints());
        assert_eq!(first.get_lookup_table(), second.get_lookup_table());

        std::fs::remove_file(nc_path).unwrap();
        std::fs::remove_file(cached).unwrap();
    }
}
