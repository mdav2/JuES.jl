using Test
using DiskTensors

tol = 1E-14
@testset "DiskTensors" begin
	@testset "Smoke" begin #checking constructors
		@testset "DiskVector" begin
			N = 1000
			x = DiskVector("/tmp/testdv1.h5",Float64,N,"w")
			y = DiskVector("/tmp/testdv2.h5",Float64,N,"w")
			@test true
		end
		@testset "DiskMatrix" begin
			N = 1000
			a = DiskMatrix("/tmp/matrix1.h5",Float64,N,N,"w")
			b = DiskMatrix("/tmp/matrix2.h5",Float64,N,N,"w")
			@test true
		end
		@testset "DiskFourTensor" begin
			DiskFourTensor("/tmp/dtens1.h5",Float64,10,10,10,10,"w")
			@test true
		end
	end
	@testset "Dot" begin
		@testset "DiskVector" begin
			N = 20000
			x = DiskVector("/tmp/testdv1.h5",Float64,N,"w")
			y = DiskVector("/tmp/testdv2.h5",Float64,N,"w")
			blockfill!(x,1.0)
			blockfill!(y,1.0)
			@test (abs(dvdot(x,y,5000) - 2E4) < tol) && (abs(dvdot(x,y) - 2E4 < tol))
			dvwrite!(x,2.0,:)
			@test abs(dvdot(x,y) - 4E4) < tol
		end
		@testset "DiskMatrix" begin
			N = 1000
			a = DiskMatrix("/tmp/matrix1.h5",Float64,N,N,"w")
			b = DiskMatrix("/tmp/matrix2.h5",Float64,N,N,"w")
			blockfill!(a,1.0)
			blockfill!(b,2.0)
			dmdot(a,b,201) - 2E6
			@test abs(dmdot(a,b,201) - 2E6) < tol
		end
		@testset "DiskFourTensor" begin
			DiskFourTensor("/tmp/dtens1.h5",Float64,10,10,10,10,"w")
			@test true
		end
	end
end
