#include "spatraster.h"
#include "string_utils.h"


SpatRaster::SpatRaster(std::string fname) {
	constructFromFile(fname);
}

SpatRaster::SpatRaster(std::vector<std::string> fname) {
	constructFromFile(fname[0]);
	SpatRaster r;
	bool success;
	for (size_t i=1; i<fname.size(); i++) {
		success = r.constructFromFile(fname[i]);
		if (success) {
			addSource(r);
			if (r.msg.has_error) {
				setError(r.msg.error);
				return;
			}
		} else {
			if (r.msg.has_error) {
				setError(r.msg.error);
			}
			return;
		}	
	}
}


void SpatRaster::setSources(std::vector<RasterSource> s) {
	source = s;
	nrow = s[0].nrow;
	ncol = s[0].ncol;
	extent = s[0].extent;
	crs = s[0].crs;
}


void SpatRaster::setSource(RasterSource s) {
	s.resize(s.nlyr);
	setSources({s});
}


SpatRaster::SpatRaster(RasterSource s) {
	setSources( {s} );
}


SpatRaster::SpatRaster() {

	RasterSource s;
	s.nrow = 10;
	s.ncol = 10;
	s.extent = SpatExtent();
	s.memory = true;
	s.filename = "";
	s.driver = "";
	s.nlyr = 1; // or 0?
	s.hasRange = { false };
	s.hasValues = false;
	s.layers.resize(1,1);
	s.datatype = "";
	s.names = {"lyr.1"};
	s.crs = "+proj=longlat +datum=WGS84";

	setSource(s);
}


SpatRaster::SpatRaster(std::vector<unsigned> rcl, std::vector<double> ext, std::string _crs) {

	RasterSource s;
	s.nrow=rcl[0];
	s.ncol=rcl[1];
	s.extent.xmin = ext[0];
	s.extent.xmax = ext[1];
	s.extent.ymin = ext[2];
	s.extent.ymax = ext[3];
	s.hasValues = false;
	s.hasRange = {false};

	s.memory = true;
	s.filename = "";
	s.driver = "";
	s.nlyr = rcl[2];
	s.layers.resize(1, 1);
	s.datatype = "";
	s.crs =_crs;
	for (unsigned i=0; i < rcl[2]; i++) { s.names.push_back("lyr." + std::to_string(i+1)) ; }

	setSource(s);
}


SpatRaster::SpatRaster(unsigned _nrow, unsigned _ncol, unsigned _nlyr, SpatExtent ext, std::string _crs) {

	RasterSource s;
	s.nrow = _nrow;
	s.ncol = _ncol;
	s.extent = ext;
	s.hasValues = false;
	s.memory = true;
	s.filename = "";
	s.driver = "";
	s.nlyr = _nlyr;
	s.hasRange = { false };
	s.layers.resize(1, 1);
	s.datatype = "";
	s.crs=_crs;
	for (unsigned i=0; i < _nlyr; i++) {	s.names.push_back("lyr." + std::to_string(i+1)) ; }
	setSource(s);
}


SpatRaster SpatRaster::geometry(long nlyrs) {
	RasterSource s;
	s.nrow = nrow;
	s.ncol = ncol;
	s.extent = extent;
	s.crs = crs;
	s.memory = true;
	nlyrs = (nlyrs < 1) ? nlyr(): nlyrs;
	s.resize(nlyrs);
	
	std::vector<std::string> nms(s.nlyr);
	for (size_t i=0; i < s.nlyr; i++) { nms[i] = "lyr" + std::to_string(i+1); }
	s.names = nms;
	SpatRaster out;
	out.setSource(s);
	return out;
}


SpatRaster SpatRaster::deepCopy() {

	SpatRaster out = *this;
	out.nrow = nrow;
	out.ncol = ncol;
	out.extent = extent;
	out.crs = crs;
	return out;
}


void SpatRaster::setCRS(std::string _crs) {
	lrtrim(_crs);
	for (size_t i = 0; i < nsrc(); i++) { source[i].crs = _crs; }
	crs = _crs;
}

std::vector<double> SpatRaster::resolution() {
	return std::vector<double> { (extent.xmax - extent.xmin) / ncol, (extent.ymax - extent.ymin) / nrow };
}

unsigned SpatRaster::nlyr() {
	unsigned x = 0;
	for (size_t i=0; i<source.size(); i++) { x += source[i].nlyr; }
	return(x);
}

std::vector<std::string> SpatRaster::filenames() {
	std::vector<std::string> x(source.size());
	for (size_t i=0; i<x.size(); i++) { x[i] = source[i].filename; }
	return(x);
}

std::vector<bool> SpatRaster::inMemory() {
	std::vector<bool> m(source.size());
	for (size_t i=0; i<m.size(); i++) { m[i] = source[i].memory; }
	return(m);
}

std::vector<bool> SpatRaster::hasRange() {
	std::vector<bool> x;
	for (size_t i=0; i<source.size(); i++) { 
		x.insert(x.end(), source[i].hasRange.begin(), source[i].hasRange.end()); 
	}
	return(x);
}

std::vector<double> SpatRaster::range_min() {
	std::vector<double> x;
	for (size_t i=0; i<source.size(); i++) { 
		x.insert(x.end(), source[i].range_min.begin(),source[i].range_min.end()); 
	}
	return(x);
}

std::vector<double> SpatRaster::range_max() {
	std::vector<double> x;
	for (size_t i=0; i<source.size(); i++) { 
		x.insert(x.end(), source[i].range_max.begin(), source[i].range_max.end()); 
	}
	return(x);
}

