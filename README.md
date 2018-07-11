# toy_MC

Toy Monte Carlo simulation for Fourier harmonics. 

Total number of events 'NumEvnts' is divided by the number of tries 'ntries' in order to generate statistical uncertainties.

To run code:

Specify settings in v1MC.C (input vn values, number of events, holes in detector acceptance, etc.). 
Then run in Root: 

'root v1MC.C+'

Output file will be stored as 'results/MH_toy.root'
