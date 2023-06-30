% Generate a test dataset for lat/lon distance

ofn = fullfile(fileparts(mfilename("fullpath")), "distance.sample.nc");

n = 20000;
tbl = table();
tbl.lat0 = wrapTo180(rand(n,1) * 360) / 2;
tbl.lat1 = wrapTo180(rand(n,1) * 360) / 2;
tbl.lon0 = wrapTo180(rand(n,1) * 360);
tbl.lon1 = wrapTo180(rand(n,1) * 360);

tbl.dist = distance(tbl.lat0, tbl.lon0, tbl.lat1, tbl.lon1, wgs84Ellipsoid("meter"));

ncid = netcdf.create(ofn, "NETCDF4");
dimid = netcdf.defDim(ncid, "i", n);

names = string(tbl.Properties.VariableNames);
varid = dictionary();

for name = names
    varid(name) = netcdf.defVar(ncid, name, "NC_DOUBLE", dimid);
    netcdf.defVarDeflate(ncid, varid(name), false, true, 9);
end % for name

netcdf.endDef(ncid)

for name = names
    netcdf.putVar(ncid, varid(name), tbl.(name));
end % for name

netcdf.close(ncid);
