cppFunction(
NumericMatrix use_surface(NumericMatrix xy_wolv, NumericMatrix xy, int n_grid, double sd_x, double sd_y, double trunc_cutoff)
{	
  int pixels = xy.nrow();
  int N_wolv = xy_wolv.nrow();
  
  double * overall_use = new double[pixels];
  double tot;
	double term1, term2, term3;
  
  for(int px = 0; px < pixels; px++)
    overall_use[px] = 1;
          
  NumericMatrix out(n_grid,3);
  
  for(int wolv = 0; wolv < N_wolv; wolv++) //Builds probability of use surface for each individual.
	{
		tot = 0.0;
		for(int px = 0; px < pixels; px++)
		{
			//Bivariate normal formula, broken up.
			term1 = pow((xy(px,0) - xy_wolv[wolv,0])/ sd_x, 2) + pow((xy(px,1) - y_wolv(wolv,1))/ sd_y, 2); 
			term2 = 2* PI * sd_x * sd_y;
			term3 = log(term2);

			use[px] = exp((term1 + 2 * term3)/(-2));
			if(trunc_cutoff>0)
			{
				if(term1 > pow(*trunc_cutoff,2)) //truncates home ranges
					use[px] = 0;
			}
			use[px] = xy(px,3) * use[px]; // rescale by snow values
			tot += use[px]; //keep track of total probability over entire surface.
		}
		for(int px = 0; px < pixels; px++){
      use[px] = use[px] / tot;    //Rescale to 1 for landscape
      out(xy(px,4),0) = 
      
			overall_use[px] = overall_use[px] * (1 - use[px] / tot); //accumulates across individuals, probability of no wolverines present.
	}
  
  delete [] overall_use;
  return out;
  }

)

 use_surface(matrix(0,3,2), matrix(0,5,4), 7, 0.4,0.5,0.3)

	double *use = new double[*pixels]; //For the cumulative probability across all individuals.
	double *overall_use = new double[*pixels];
	const double PI = 3.141593;


	for(int px = 0; px < *pixels; px++)
		overall_use[px] = 1;

	for(int wolv = 0; wolv < *N_wolv; wolv++) //Builds probability of use surface for each individual.
	{
		tot = 0.0;
		for(int px = 0; px < *pixels; px++)
		{
			//Bivariate normal formula, broken up.
			term1 = pow((x[px] - x_wolv[wolv])/ *sd_x, 2) + pow((y[px] - y_wolv[wolv])/ *sd_y, 2); 
			term2 = 2* PI * *sd_x * *sd_y;
			term3 = log(term2);

			use[px] = exp((term1 + 2 * term3)/(-2));
			if(*trunc_cutoff>0)
			{
				if(term1 > pow(*trunc_cutoff,2)) //truncates home ranges
					use[px] = 0;
			}
			use[px] = snow[px] * use[px]; // rescale by snow values
			tot += use[px]; //keep track of total probability over entire surface.
		}
		for(int px = 0; px < *pixels; px++)
			overall_use[px] = overall_use[px] * (1 - use[px] / tot); //accumulates across individuals, probability of no wolverines present.
	}

	for(int px = 0; px < *pixels; px++)
		snow[px] = 1 - overall_use[px]; // 1- probability of no wolverines = probability of at least one wolverine.

	delete[] overall_use;
	delete[] use;
  
}
)