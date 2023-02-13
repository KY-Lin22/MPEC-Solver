function MPEC = checkInput(self, MPEC)
%MPEC = checkInput(self, MPEC)
%   check and polish the input of NIPMPEC_CasADi: MPEC struct
      
% check x
if isempty(MPEC.x)
    error('please specify optimal variable MPEC.x')
else
    if size(MPEC.x, 2) ~= 1
        error('MPEC.x should be a column vector')
    end
end

% check p
if isempty(MPEC.p)
    MPEC.p = casadi.SX.sym('p', 0, 1);
else
    if size(MPEC.p, 2) ~= 1
        error('MPEC.p should be a column vector')
    end
end

% check L
if isempty(MPEC.L)
    error('please specify cost function MPEC.L')
else
    if ~all(size(MPEC.L) == [1, 1])
        error('MPEC.L should be a scale function')
    end
end

% check G
if isempty(MPEC.G)
    MPEC.G = casadi.SX.sym('G', 0, 1);
else
    if size(MPEC.G, 2) ~= 1
        error('MPEC.G should be a column function')
    end
end

% check C
if isempty(MPEC.C)
    MPEC.C = casadi.SX.sym('C', 0, 1);
else
    if size(MPEC.C, 2) ~= 1
        error('MPEC.C should be a column function')
    end
end

% check K
if isempty(MPEC.K)
    if size(MPEC.p, 1) == 0
        MPEC.K = casadi.SX.sym('K', 0, 1);
    else
        error('please specify MPEC.K')
    end
else
    if ~all(size(MPEC.K) == size(MPEC.p))
        error('MPEC.K should have the same dimension as MPEC.p')
    end
end

% check l
if isempty(MPEC.l)
    if size(MPEC.p, 1) == 0
        MPEC.l = zeros(0, 1);
    else
        error('please specify MPEC.l')
    end   
else   
    MPEC.l = double(MPEC.l);
    if ~all(size(MPEC.l) == size(MPEC.p))
        error('MPEC.l should have the same dimension as MPEC.p')
    end
end

% check u
if isempty(MPEC.u)
    if size(MPEC.p, 1) == 0
        MPEC.u = zeros(0, 1);
    else
        error('please specify MPEC.u')
    end   
else    
    MPEC.u = double(MPEC.u);
    if ~all(size(MPEC.u) == size(MPEC.p))
        error('MPEC.u should have the same dimension as MPEC.p')
    end
end

% check l <= u
if ~all(MPEC.l <= MPEC.u)
    error('MPEC.l should less than or equal to MPEC.u')
end

end

