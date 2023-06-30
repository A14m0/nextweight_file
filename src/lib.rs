use std::str::FromStr;
use std::{path::PathBuf, io::Write, mem::size_of, io::Read};



use std::collections::HashMap;

use netcdf::AttrValue;

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
    per_variable_attrs: HashMap<String, Vec<(String, String)>>,
    polyids: Vec<String>
}

#[derive(Debug)]
#[repr(C)]
pub struct PolyidEntry {
    pub data: Vec<(u32, u32, f32, f32, f32)>
}

impl NextWeightFile {
    /// opens a NetCDF weight file and converts it to
    pub fn from_weight_file(path: PathBuf) -> Result<Self, String> {
        // open the weight file
        let weight_netcdf = netcdf::open(path).unwrap();
        let mut json_data = JsonData::new();

        // now we get all of the attributes...
        for attr in weight_netcdf.attributes() {
            let attr_value = match attr.value().unwrap() {
                AttrValue::Str(a) => a,
                AttrValue::Strs(a) => a[0].clone(),
                _ => return Err(format!("Unexpected attribute type for {}", attr.name()))
            };
            // add it to our list of global attributes
            json_data.add_global_attr(attr.name().to_string(), attr_value);
        }

        // ... and add all of the variable attributes ...
        for var in weight_netcdf.variables() {
            let var_name = var.name();
            json_data.add_variable(&var_name);
            for attr in var.attributes() {
                if attr.name() != "_FillValue" {
                    let attr_value = match attr.value().unwrap() {
                        AttrValue::Str(a) => a,
                        AttrValue::Strs(a) => a[0].clone(),
                        _ => return Err(format!("Unexpected attribute type for {}", attr.name()))
                    };
                    // add it to our list of global attributes
                    json_data.add_variable_attr(&var_name, attr.name().to_string(), attr_value);
                }
            }
        }

        // now that we have gotten our attributes all squared away, lets start
        // looking at data. First things first, lets store those polyids
        let polyid_var = weight_netcdf.variable("polyid").unwrap();
        for polyid in 0..polyid_var.len() {
            json_data.add_polyid(polyid_var.string_value(polyid).unwrap());
        };

        // next lets start processing those weights
        let regridweights = weight_netcdf.variable("regridweights").unwrap();
        let latvar = weight_netcdf.variable("lat").unwrap();
        let lonvar = weight_netcdf.variable("lon").unwrap();
        let lat_vals = latvar.values::<f32,_>(..).unwrap();
        let lon_vals = lonvar.values::<f32,_>(..).unwrap();
        let lat_len = weight_netcdf.dimension("lat").unwrap().len() as u64;
        let lon_len = weight_netcdf.dimension("lon").unwrap().len() as u64;
        let fill = regridweights.fill_value::<f32>().unwrap().unwrap();
        let mut polyid_gridpoints: Vec<PolyidEntry> = Vec::new();

        // for every polyid...
        for polyid in 0..polyid_var.len() {
            // ... create a new entry into our lookup vector...
            let mut curr_polyid = PolyidEntry::new();
            let data = regridweights.values_arr::<f32,_>((polyid,..,..)).unwrap();
            let dat_slice = data.as_slice().unwrap();
            // ... for every data value...
            for lat_idx in 0..lat_len as usize {
                for lon_idx in 0..lon_len as usize {
                    let data_value = dat_slice[lat_idx * lon_len as usize+  lon_idx];
                    // ...if it isnt a fill value...
                    if data_value != fill {
                        // ... then calculate the lat lon and save the weight
                        curr_polyid.add_point(lat_idx as u32, lon_idx as u32, lat_vals[lat_idx], lon_vals[lon_idx], data_value);
                    }
                }
            }

            // now push the polyid entry to our lookup vector
            polyid_gridpoints.push(curr_polyid);
        }

        // and finally lets build our lookup table
        let mut lookup_table: Vec<(u64, u64)> = Vec::new();
        let mut running_total: u64 = 0;
        for entry in polyid_gridpoints.iter() {
            let entry_size = entry.data.len()as u64;
            lookup_table.push((running_total, entry.data.len() as u64));
            running_total += entry_size;
        }

        // now we are done, so return ourselves
        Ok(Self {
            json_data,
            lat_len,
            lon_len,
            polyid_gridpoints,
            lookup_table
        })
    }

    /// create new structure from .NWT file
    pub fn from_nwt(path: PathBuf) -> Result<Self, String> {
        // open the file
        let mut input_file = std::fs::File::open(path).unwrap();
        let mut data: Vec<u8> = Vec::new();
        input_file.read_to_end(&mut data).unwrap();

        let mut co: usize = 0;
        // first check for magic
        if &data[co..co+4] != b"NEWT" {
            return Err("Invalid file format".to_string());
        }
        co += 4;

        // now we read all the crap we need
        let mut u64_buff = [0u8; size_of::<u64>()];
        // json_strlen
        u64_buff.copy_from_slice(&data[co..co+8]);
        let json_len = u64::from_le_bytes(u64_buff);
        co += 8;
        // number of polyids
        u64_buff.copy_from_slice(&data[co..co+8]);
        let num_polyids = u64::from_le_bytes(u64_buff);
        co += 8;
        // latitude length
        u64_buff.copy_from_slice(&data[co..co+8]);
        let lat_len = u64::from_le_bytes(u64_buff);
        co += 8;
        // longitude length
        u64_buff.copy_from_slice(&data[co..co+8]);
        let lon_len = u64::from_le_bytes(u64_buff);
        co += 8;
        // json string offset
        u64_buff.copy_from_slice(&data[co..co+8]);
        let json_offset = u64::from_le_bytes(u64_buff);
        co += 8;
        // lookup offset
        u64_buff.copy_from_slice(&data[co..co+8]);
        let lookup_offset = u64::from_le_bytes(u64_buff);
        
        // json data
        let mut string_buffer: Vec<u8> = vec![0u8;json_len as usize];
        string_buffer.copy_from_slice(&data[(json_offset as usize )..(json_offset+json_len) as usize]);
        let json_dat_str = String::from_utf8(string_buffer).unwrap();
        let json_data = serde_json::from_str(&json_dat_str).unwrap();
        
        // now we get the lookup table information
        let mut lookup_table_bytes: Vec<u8> = vec![0u8;num_polyids as usize * size_of::<(u64,u64)>()]; //Vec::with_capacity(num_polyids as usize * size_of::<u64>());
        lookup_table_bytes.copy_from_slice(
            &data[
                (lookup_offset as usize) ..
                (lookup_offset as usize + num_polyids as usize *size_of::<(u64,u64)>())
                ]);

        let mut lookup_table: Vec<(u64,u64)> = Vec::new();
        for i in (0..lookup_table_bytes.len()).step_by(16) {
            u64_buff.copy_from_slice(&lookup_table_bytes[i..i+8]);
            let tmp_u64 = u64::from_le_bytes(u64_buff);
            u64_buff.copy_from_slice(&lookup_table_bytes[i+8..i+16]);
            let count_u64 = u64::from_le_bytes(u64_buff);
            lookup_table.push((tmp_u64, count_u64));
        }
        
        // and finally now that we have that, we pull all of our weight values
        let mut polyid_gridpoints: Vec<PolyidEntry> = Vec::new();
        co = lookup_offset as usize + num_polyids as usize * size_of::<(u64,u64)>();
        for i in 0..num_polyids {
            // read in the number of grid coordinates we are to expect
            let num_coords = lookup_table[i as usize].1;
            let mut buff_32 = [0u8; size_of::<f32>()];
            let mut curr_polyid = PolyidEntry::new();
            for _ in 0..num_coords {
                // read lat index
                buff_32.copy_from_slice(&data[co..co+4]);
                co += 4;
                let lat_idx = u32::from_le_bytes(buff_32);
                // read lon index
                buff_32.copy_from_slice(&data[co..co+4]);
                co += 4;
                let lon_idx = u32::from_le_bytes(buff_32);
                // read actual latitude
                buff_32.copy_from_slice(&data[co..co+4]);
                co += 4;
                let lat = f32::from_le_bytes(buff_32);
                // read actual longitude
                buff_32.copy_from_slice(&data[co..co+4]);
                co += 4;
                let lon = f32::from_le_bytes(buff_32);
                // read weight
                buff_32.copy_from_slice(&data[co..co+4]);
                co += 4;
                let weight = f32::from_le_bytes(buff_32);

                // and add it to our list
                curr_polyid.add_point(lat_idx, lon_idx, lat, lon, weight);
            }

            // add the polyid to our polyid gridpoitns
            polyid_gridpoints.push(curr_polyid);
        }


        // now that we have everything, lets return stuff

        Ok(Self { json_data, lat_len, lon_len, polyid_gridpoints, lookup_table })


    }

    /// Generically opens a weight file. If it is a NetCDF file, it is converted 
    /// to the NWT format. Otherwise it is opened as standard
    pub fn open(path: PathBuf) -> Result<Self, String> {
        let mut data = [0u8; 4];
        // scope brackets here to make sure `input_file` is closed before opening
        {
            let mut input_file = std::fs::File::open(&path).unwrap();
            input_file.read(&mut data).unwrap();
        }

        // first check for magic
        if &data[0..4] == b"NEWT" {
            Self::from_nwt(path)
        } else {
            let t_path = path.clone();
            let name = t_path.to_str().unwrap();
            let a = Self::from_weight_file(path)?;
            let new_path = PathBuf::from_str(&format!("{}.nwt", name)[..]).unwrap();
            println!("[libNextWeightFile] Serializing new weight file to {}. Use this next time to avoid precomputation step", new_path.display());
            a.serialize_to_file(Some(new_path.to_str().unwrap().to_string()))?;
            Ok(a)
        }
    }

    /// serializes the new weight file to disk
    pub fn serialize_to_file(&self, filename: Option<String>) -> Result<(), String> {
        // first determine our filename. Default is "test.nwt"
        let fname = match filename {
            Some(a) => a,
            None => "test.nwt".to_string()
        };
        // then lets create/open our file
        let mut output_file = match std::fs::File::create(fname) {
            Ok(a) => a,
            Err(e) => return Err(format!("Failed to open file for serialization: {}",e))
        };

        // first we write some of the important things we need in the header
        let serialized_dat = serde_json::to_string(&self.json_data).unwrap();
        // magic bytes
        output_file.write(b"NEWT").unwrap();
        // u64: length of json string
        output_file.write(&(serialized_dat.len() as u64).to_le_bytes()).unwrap();
        // u64: number of polyids
        output_file.write(&(self.json_data.polyids.len() as u64).to_le_bytes()).unwrap();
        // u64: latitude length
        output_file.write(&self.lat_len.to_le_bytes()).unwrap();
        // u64: longitude length
        output_file.write(&self.lon_len.to_le_bytes()).unwrap();
        // beginning of json attributes string
        let json_offset = size_of::<u64>() * 6 + 4;
        output_file.write(&json_offset.to_le_bytes()).unwrap();
        // beginning of lookup vector
        let lookup_offset = json_offset + serialized_dat.len();
        output_file.write(&lookup_offset.to_le_bytes()).unwrap();
        // the actual json data
        write!(output_file, "{}", serialized_dat).unwrap();

        // next we build our lookup table
        for v in self.lookup_table.iter() {
            output_file.write(&v.0.to_le_bytes()).unwrap();
            output_file.write(&v.1.to_le_bytes()).unwrap();
        }

        // and finally we can now serialize all data
        for d in self.polyid_gridpoints.iter() {
            // and then the values
            for v in d.data.iter() {
                output_file.write(&v.0.to_le_bytes()).unwrap();
                output_file.write(&v.1.to_le_bytes()).unwrap();
                output_file.write(&v.2.to_le_bytes()).unwrap();
                output_file.write(&v.3.to_le_bytes()).unwrap();
                output_file.write(&v.4.to_le_bytes()).unwrap();
            }
        }

        Ok(())

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
    pub fn get_raw_gridpoints(&self) -> Vec<(u32, u32, f32, f32, f32)> {
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


impl JsonData {
    /// creates a new instance of `JsonData`
    pub fn new() -> Self {
        Self {
            global_attrs: Vec::new(),
            per_variable_attrs: HashMap::new(),
            polyids: Vec::new()
        }
    }

    /// adds a global attribute to the structure
    pub fn add_global_attr(&mut self, key: String, value: String) {
        self.global_attrs.push((key, value));
    }

    /// adds a new variable to the structure
    pub fn add_variable(&mut self, variable_name: &String) {
        self.per_variable_attrs.insert(variable_name.clone(), Vec::new());
    }

    /// adds a new attribute for the associated variable. If the variable has
    /// not yet been added, it is added
    pub fn add_variable_attr(&mut self, var_name: &String, key: String, value: String) {
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
    pub fn get_global_attr(&self, name: &String) -> Result<String, String> {
        for v in self.global_attrs.iter() {
            if v.0 == *name { return Ok(v.1.clone()) }
        }
        Err(format!("Global attribute {} not found", name))
    }

    /// Retrieves all global attributes
    pub fn get_global_attrs(&self) -> &Vec<(String, String)> {
        &self.global_attrs
    }

    /// Retrieves a given variable's attribute of a provided name
    pub fn get_var_attr(&self, variable_name: &String, attr_name: &String) -> Result<String, String> {
        match self.per_variable_attrs.get(variable_name) {
            Some(a) => {
                for v in a.iter() {
                    if v.0 == *attr_name { return Ok(v.1.clone()) }
                }
                Err(format!("Global attribute {} not found", attr_name)) 
            },
            None => return Err(format!("Variable {} not found in the weight file", variable_name))
        }
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



#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn it_works() {
        // lets test this
        let test_path = PathBuf::from_str("test_cases/cckp_aggregation_1x1").unwrap();
        let new_path = PathBuf::from_str("test.nwt").unwrap();
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
}
