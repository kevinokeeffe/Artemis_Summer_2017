classdef PropFund < handle
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties
		
		%	Fundamental wavelength [nm] (can be non-scalar)
		l0 double {mustBePositive, mustBeFinite, mustBeNonempty, ...
			mustBeReal} = 780;
		
		%	Focal length [mm] (must be scalar)
		f(1,1) double {mustBePositive, mustBeFinite, mustBeNonempty, ...
			mustBeReal} = 300;
		
	end %	properties
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (SetAccess=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (SetAccess=private)
		H = hankel_matrix(0, 1, 2^8);					%	Hankel matrix
	end %	properties

	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (Dependent, SetAccess=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (Dependent)
		r_FF;			%	Radial co-ordinates in far field/focus [mm]
		kr_FF;			%	Radial wavevector in far field/focus [rad/mm]
		kz_FF;			%	Longitudinal wavevector in far field/focus [rad/mm]
		r_NF;			%	Radial co-ordinates in near field [mm]
		kr_NF;			%	Radial wavevector in near field [rad/mm]
		kz_NF;			%	Longitudinal wavevector in near field/focus [rad/mm]
		k0;				%	Fundamental wavevector [rad/mm]
	end
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	--------------------------------------------------------------------
		
		function obj = PropFund(varargin)
			
			%	Set up Hankel matrix
			obj.SetHankelMatrix(varargin);			
		end %	function obj = HHG_QPI_Simple(varargin)
	
		%	--------------------------------------------------------------------
		
		function SetHankelMatrix(obj, varargin)
			Nr = FindArg(varargin, 'Nr', obj.H.Nr);
			Rmax = FindArg(varargin, 'RMax', obj.H.R);
			
			%	Only generate hankel matrix if necessary
			if ~isequal(Nr, obj.H.Nr) || ~isequal(Rmax, obj.H.R)
				obj.H = hankel_matrix(0, Rmax, Nr);
			end
			
		end %	function SetHankelMatrix(obj, varargin)
		
		%	--------------------------------------------------------------------
		
		%	Calculates the far-field (i.e. focused) electric-field distribution
		%	given the near-field input
		function E_FF = CalculateFarField(obj, E_NF)
			%	Need to scale field before Hankel transform
			U = E_NF ./ obj.H.JR;
			
			%	Perform Hankel transform to calculate far-field
			U = reshape(obj.H.T * U(:, :), size(U));
			
			%	Rescale field back to "physical" space
			E_FF = obj.H.JV .* U;
		end
		
		%	--------------------------------------------------------------------
		
		function Ez = PropagateFF(obj, E_FF, z)
			Ez = obj.Propagate(E_FF, z, obj.kz_FF);
		end
		
		function Ez = PropagateNF(obj, E_FF, z)
			Ez = obj.Propagate(E_FF, z, obj.kz_NF);
		end
		
		%	--------------------------------------------------------------------

	end %	methods
	

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (set/get)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

		%	--------------------------------------------------------------------
		
		%	Get fundamental wavevector [rad/mm]
		function val = get.k0(obj)
			val = 2e6*pi./obj.l0;
		end
		
		%	Get radial co-ordinates at focus [mm]
		function val = get.r_FF(obj)
			val = obj.H.kr .* (obj.f ./ obj.k0);
		end
		
		%	Get radial wavevector at focus [rad/mm]
		function val = get.kr_FF(obj)
			val = obj.H.r .* (obj.k0 ./ obj.f);
		end
		
		%	Get radial co-ordinates at focus [mm]
		function val = get.r_NF(obj)
			val = obj.H.r;
		end
		
		%	Get radial wavevector at focus [rad/mm]
		function val = get.kr_NF(obj)
			val = obj.H.kr;
		end
		
		%	Get longitudinal wavevector in near field [rad/mm]
		function val = get.kz_NF(obj)
			val = sqrt((obj.k0.^2 - obj.H.kr.^2) .* (obj.k0>obj.H.kr)) - obj.k0;
		end
		
		%	Get longitudinal wavevector at focus [rad/mm]
		function val = get.kz_FF(obj)
			val = sqrt((obj.k0.^2 - obj.kr_FF.^2) .* (obj.k0>obj.kr_FF)) - obj.k0;
		end
		
		%	--------------------------------------------------------------------
		
	end %	methods (set/get)
		
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (Access=protected)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Access=protected)
		
		%	--------------------------------------------------------------------

		function Ez = Propagate(obj, E, z, kz)
			NE = size(E);
			Nz = length(z);
			N = [NE Nz];
			N(N<=1) = [];
			if isscalar(N)
				N = [N 1];
			end
			
			%	Need to scale field before Hankel transform
			U = E ./ obj.H.JR;
			
			%	Perform Hankel transform
			U = obj.H.T * U(:, :);
			
			%	Apply propagation phase factor
			U = U .* exp(1i*kz.*reshape(z, 1, Nz));
			
			%	Perform Hankel transform
			U = reshape(obj.H.T * U, N);
			
			%	Rescale field
			Ez = U .* obj.H.JV;
			
		end %	function Ez = Propagate(...)
		
		%	--------------------------------------------------------------------

	end %	methods (Access=protected)
	
	
end	%	classdef
