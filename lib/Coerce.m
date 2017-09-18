%	Coerce values into defined range
%
%	[out, coerced] = Coerce(in, [mn], [mx]
%
%	If mn defined:
%		out(in<mn) = mn
%
%	If mx defined:
%		out(in>mx) = mx
%
%	Author(s)
%		Dr Adam S Wyatt (adam.wyatt@stfc.ac.uk)
function [out, coerced] = Coerce(in, mn, mx)

if nargin>=3 && ~isempty(mn) && ~isempty(mx)
	if mx<mn
		[mx, mn] = Swap(mn, mx);
	end
end

out = in;
coerced = false(size(out));

if nargin>=2 && ~isempty(mn)
	ind = (in<mn);
	if isscalar(mn)
		out(ind) = mn;
	else
		out(ind) = mn(ind);
	end
	coerced(ind) = true;
end

if nargin>=3 && ~isempty(mx)
	ind = (in>mx);
	if isscalar(mx)
		out(ind) = mx;
	else
		out(ind) = mx(ind);
	end
	coerced(ind) = true;
end

end
