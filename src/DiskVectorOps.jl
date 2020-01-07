function dvdot(dvec1::DiskVector,dvec2::DiskVector)
	"computes dot product of two DiskVectors - elementwise, very slow!"
	if (dvec1.size != dvec2.size)
		return false
	end
	tsum = 0.0
	for i in 1:1:dvec1.size
		tsum += dvec1[i]*dvec2[i]
	end
	return tsum
end
function dvdot(dvec1::DiskVector,dvec2::DiskVector,buffsize)
	"computes dot product of two DiskVectors - buffered, use this"
	chunks = cld(dvec1.size,buffsize)
	tsum = 0.0
	for chunk in 1:1:(chunks)
        tsum += dot(dvec1[(chunk-1)*buffsize+1:(chunk)*buffsize],
                    dvec2[(chunk-1)*buffsize+1:(chunk)*buffsize])
				
	end
	return tsum
end
