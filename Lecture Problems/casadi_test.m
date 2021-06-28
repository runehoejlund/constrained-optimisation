addpath('/Users/rune.hoejlund/Development/DTU/Kandidat/casadi-osx-matlabR2015a-v3.5.5')
import casadi.*
x = MX.sym('x')
disp(jacobian(sin(x),x))
