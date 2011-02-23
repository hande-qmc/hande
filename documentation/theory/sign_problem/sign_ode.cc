#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cctype>
#include <cmath>

using namespace std;

void strip_comments(string &line);
string decommented_file_as_string(const char *filename);

class Rhs {
  double a;
  double b;
public:
  Rhs(double aa, double bb) : a(aa), b(bb) {};
  Rhs(const Rhs& rhsin) : a(rhsin.a), b(rhsin.b) {};
  double operator()(double t, double x) const {
    return x - a*x*x + exp(b*t);
  }
};

double rk4increment(const Rhs& rhs, double t, 
		    double x, double dt);

int main()
{
  // This program solves the ODE
  //     dx/dt = x - ax^2 + e^{bt}
  // with the initial condition x(0) = 0.

  double a, b, dt, nsteps;
  double t = 0.0;
  double x = 0.0;
  
  string decommented_file = decommented_file_as_string("sign_ode.dat");
  istringstream instream(decommented_file);
  instream >> a >> b >> dt >> nsteps;
  Rhs rhs(a, b);

   cout << t << '\t' << x << '\n';
  for (int istep=0; istep < nsteps; ++istep) {
    x += rk4increment(rhs, t, x, dt);
    t += dt;
    cout << t << '\t' << x << '\n';
  }

  return 0;
}

string decommented_file_as_string(const char *filename)
{
  ifstream infile(filename);
  string line;
  string decommented_file;
  while ( getline(infile, line, '\n') ) {
    strip_comments(line);
    if (line.size() != 0) {
      decommented_file += line += '\n';
    }
  }
  return decommented_file;
}

void strip_comments(string &line)
{
  // Strip everything from the first '#' character
  // to the end of the input string.
  string::size_type pos = 0;
  if ( (pos = line.find('#')) != string::npos ){
    line.erase(pos, string::npos);
  }
  // Strip trailing spaces, if any.
  pos = line.size() - 1;
  while ( isspace(line[pos]) ) --pos;
  line.erase(pos+1);
}

double rk4increment(const Rhs& rhs, double t, 
		    double x, double dt)
{
  // Use 4th order Runge-Kutta to move forward by dt.
  double k1 = rhs(t, x);
  double k2 = rhs(t + 0.5*dt, x + 0.5*dt*k1);
  double k3 = rhs(t + 0.5*dt, x + 0.5*dt*k2);
  double k4 = rhs(t + dt, x + dt*k3);
  return dt * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
}
