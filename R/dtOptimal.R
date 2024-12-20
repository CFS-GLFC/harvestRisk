## Given a numeric vector, return most efficient data type to write using raster::writeRaster()

dtOptimal <- function(x) {
  
  rval <- range(x)
  
  if(all(rval %% 1 == 0)) {
    
    if(rval[1]==0 & rval[2]==1) return('LOG1S')
    if(rval[1] >= 0 & rval[2] <= 255) return('INT1U')
    if(rval[1] >= 0 & rval[2] <= 65534) return('INT2U')
    if(rval[1] >= 0 & rval[2] <= 4294967296) return('INT4U')
    if(rval[1] >= -127 & rval[2] <= 127) return('INT1S')
    if(rval[1] >= -32767 & rval[2] <= 32767) return('INT2S')
    if(rval[1] >= -2147483647 & rval[2] <= 2147483647) return('INT4S')
    
  } else {
    
    if(rval[1] >= -3.4e38 & rval[2] <= 3.4e38) return('FLT4S')
    if(rval[1] >= -1.7e308 & rval[2] <= 1.7e308) return('FLT8S')
    
  }
  
}
