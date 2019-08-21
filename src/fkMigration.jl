module fkMigration
	using MAT
	using FFTW
	
	export demo, loadDataset, calibrate!, downsampleAndCrop, reconstruct
	
	
	const c = 299792458  # speed of light (m/s)
	
	
	demo(dataset = "teaser", resXY = 512, cropT = 1024) = reconstruct(downsampleAndCrop(calibrate!(loadDataset(dataset)...), resXY, cropT))
	
	
	function loadDataset(dataset, filename = "meas_180min.mat")
		data = matread(joinpath(dataset, filename))
		tau = data["meas"]
		
		data = matread(joinpath(dataset, "tof.mat"))
		calib = data["tofgrid"]
		
		return (tau, calib)
	end
	
	
	function calibrate!(tau, calib; wallSize = 2, resT = 32e-12)
		# calibrate so that t=0 is when light reaches the scan surface
		for x in 1:size(tau, 1), y in 1:size(tau, 2)
			delay = floor(Int, calib[x, y]*1e-12 / resT)
			tau[x, y, :] .= circshift(tau[x, y, :], -delay)
		end
		return tau
	end
	
	
	function downsampleAndCrop(tau, resXY, cropT)
		function downsample(A)
			Isize = (size(A, 1)+1)÷2
			Jsize = (size(A, 2)+1)÷2
			for i in 1:Isize, j in 1:Jsize
				I = 2*(i-1) + 1
				J = 2*(j-1) + 1
				A[i, j, :] .= 0.25 .* (
					A[min(I,   size(A, 1)), min(J,   size(A, 2)), :] .+
					A[min(I+1, size(A, 1)), min(J,   size(A, 2)), :] .+
					A[min(I+1, size(A, 1)), min(J+1, size(A, 2)), :] .+
					A[min(I,   size(A, 1)), min(J+1, size(A, 2)), :]
				)
			end
			return A[1:Isize, 1:Jsize, :]
		end
		
		tau = tau[:, :, 1:min(cropT, size(tau, 3))]
		while size(tau, 1) > resXY || size(tau, 2) > resXY
			tau = downsample(tau)
		end
		return tau
	end
	
	
	function reconstruct(tau; wallSize = 2, resT = 32e-12)
		(X, Y, T) = size(tau)
		Psi = zeros(eltype(tau), map(s -> 2*s, size(tau)))
		for t in 1:T
			Psi[1:X, 1:Y, t] .= sqrt.(tau[1:X, 1:Y, t])
		end
		
		PhiBar = fftshift(fft(Psi))
		
		Phi = zero(PhiBar)
		for K_x in 1:size(Phi, 1), K_y in 1:size(Phi, 2), K_z in 1:size(Phi, 3)
			Phi[K_x, K_y, K_z] = stoltInterp(PhiBar, K_x, K_y, K_z, physWidth = wallSize, physHeight = wallSize, resT = resT)
		end
		
		Psi = ifft(ifftshift(Phi))
		Psi  = abs2.(Psi[1:X, 1:Y, 1:T])
		
		return Psi
	end
	
	
	function stoltLinearInterpZ(A, x, y, z)
		zLow  = floor(Int, z)
		zHigh = ceil(Int, z)
		low   = 1 <= zLow  <= size(A, 3) ? A[x, y, zLow] : 0
		high  = 1 <= zHigh <= size(A, 3) ? A[x, y, zHigh] : 0
		return low*(zHigh-z) + high*(z-zLow)
	end
	
	
	function stoltInterp(A, x, y, z; physWidth, physHeight, resT)
		# normalize indices and convert i and j to units of k
		i = 2*(x-1)/size(A, 1) - 1
		j = 2*(y-1)/size(A, 2) - 1
		k = 2*(z-1)/size(A, 3) - 1
		
		i *= (size(A, 1)/2*c*resT) / (4*(physWidth/2))
		j *= (size(A, 2)/2*c*resT) / (4*(physHeight/2))
		
		l = sqrt(i^2 + j^2 + k^2)
		f = (l + 1)*size(A, 3)/2 + 1
		S = max(c/2 * k/max(l, eps()), 0)
		return k > 0 ? S*stoltLinearInterpZ(A, x, y, f) : 0
	end
end
