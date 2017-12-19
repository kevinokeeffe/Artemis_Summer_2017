classdef HHG_QPI_Simple < handle
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties
		%	Harmonic order (can be non-scalar)
		q double {mustBePositive, mustBeFinite, mustBeNonempty, mustBeReal, ...
			mustBeInteger, HHG_QPI_Simple.mustBeOdd} = 21;
		
		%	Effective nonlinearity (can be non-scalar)
		q_eff double {mustBePositive, mustBeFinite, mustBeNonempty, ...
			mustBeReal} = 9;
		
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
		r;				%	Radial co-ordinates in focus [mm]
		kr;				%	Radial wavevector in focus [rad/mm]
		k0;				%	Fundamental wavevector [rad/mm]
		kq;				%	Harmonic wavevector [rad/mm]
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods
		
		%	--------------------------------------------------------------------
		
		function obj = HHG_QPI_Simple(varargin)
			
			%	Set up Hankel matrix
			obj.SetHankelMatrix(varargin);
			
			%	Set up propagation parameters
			obj.q = FindArg(varargin, 'HarmonicOrder', obj.q);
			obj.q_eff = FindArg(varargin, 'EffectiveNonlinearity', obj.q_eff);
			
		end %	function obj = HHG_QPI_Simple(varargin)
	
		%	--------------------------------------------------------------------
		
		function SetHankelMatrix(obj, varargin)
			Nr = FindArg(varargin, 'Nr', obj.H.Nr);
			Rmax = FindArg(varargin, 'RMax', obj.H.R);
			
			%	Only generate hankel matrix if necessary
			if ~isequal(Nr, obj.H.Nr) && ~isequal(Rmax, obj.H.R)
				obj.H = hankel_matrix(0, Rmax, Nr);
			end
			
		end %	function SetHankelMatrix(obj, varargin)
		
		%	--------------------------------------------------------------------
		%	--------------------------------------------------------------------
		%	--------------------------------------------------------------------

	end %	methods
	

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (set/get)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

		%	--------------------------------------------------------------------
		
		%	Get fundamental wavevector
		function val = get.k0(obj)
			val = 2*pi./obj.l0;
		end
		
		%	Get harmonic wavevector
		function val = get.kq(obj)
			val = obj.q.*obj.k0;
		end
		
		%	Get radial co-ordinates at focus
		function val = get.r(obj)
			val = obj.H.kr .* (obj.f ./ obj.k0);
		end
		
		%	Get radial wavevector at focus
		function val = get.kr(obj)
			val = obj.H.r .* (obj.k0 ./ obj.f);
		end
		
		%	--------------------------------------------------------------------
		
	end
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (Static)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods (Static)
		
		%	--------------------------------------------------------------------
		
		function mustBeOdd(val)
			if ~all(bitget(val, 1))
				error('Value must be odd');
			end		
		end
		
		%	--------------------------------------------------------------------

	end
	
end	%	classdef
