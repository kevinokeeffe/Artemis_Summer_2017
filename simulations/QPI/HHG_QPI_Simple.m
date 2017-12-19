classdef HHG_QPI_Simple < PropFund

	
	
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
	end %	properties
	
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	P R O P E R T I E S (Dependent, SetAccess=private)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties (Dependent)
		kq;				%	Harmonic wavevector [rad/mm]
	end %	properties (Dependent)
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

		%	--------------------------------------------------------------------

		function obj = HHG_QPI_Simple(varargin)
			
			%	Call superclass constructor
			obj@PropFund(varargin);
			
			
			%	Set up propagation parameters
			obj.q = FindArg(varargin, 'HarmonicOrder', obj.q);
			obj.q_eff = FindArg(varargin, 'EffectiveNonlinearity', obj.q_eff);
			
		end %	function obj = HHG_QPI_Simple(varargin)

		%	--------------------------------------------------------------------
	
	end %	methods
		

	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	M E T H O D S (set/get)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods

		%	--------------------------------------------------------------------
		
		%	Get harmonic wavevector
		function val = get.kq(obj)
			val = obj.q.*obj.k0;
		end
		
		%	--------------------------------------------------------------------
		
	end %	methods (set/get)


	
	
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

	end %	methods (Static)
	
end %	classdef
