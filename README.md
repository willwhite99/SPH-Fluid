# SPH-Fluid
My final year project at university was a SPH fluid solver made in C++, focus was given to the implementation to make it stable and real-time. Of which the solution is parallelised, uses SIMD and has a static grid for broad phase collision.

It does however have a naive implementation of SIMD that does not benefit the application that much due to using a AoS, which is not very compaitable with SIMD, however this can be changed. Also the static grid could be supplemented with a hash tree for better performance.

# References

Chentanez, N. and Müller. M. (2011) ‘Real-time Eulerian water simulation using a restricted tall cell grid.’ Proceedings of ACM SIGGRAPH 2011, Volume 30, Issue 4, Article 82

Friedrich, R. and Gurevich, S. (2010) Numerical Methods for Complex Systems. Available at: http://pauli.uni-muenster.de/tp/fileadmin/lehre/NumMethoden/WS1011/script1011.pdf (Accessed: 13 March 2015)

Gresho, P. (1990) ‘Some Current CFD Issues Relevant to the Incompressible Navier-Stokes Equations’, Journal of Computer Methods in Applied Mechanics and Engineering, Volume 87, pp. 201-252

McAdams, A., Sifakis, E. and Teran, J. (2010) ‘A parallel multigrid Poisson solver for fluids simulation on large grids’, Proceedings of the 2010 ACM SIGGRAPH/Eurographics Symposium on Computer Animation, Madrid, Spain, 2-4 July 
