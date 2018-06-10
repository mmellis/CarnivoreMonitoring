cppFunction(
NumericVector sample_ind(NumericMatrix xy, int N, double buffer, IntegerVector use, IntegerVector new_order)
{
  int maxid = use.size();
  int longlat;
  double (*rdist)(double, double, double, double);       
  longlat=use[new_order[0]];  // Checks whether to use rdist_earth or rdist to calculate distances
  if(longlat == 1)
    rdist = &rdist_earth;
  else
    rdist = &rdist_km;  
    
	use[new_order[0]] = 1; // Automatically include the first point in the sample = 1st wolverine
	int ct = 1; // Number of individuals included so far
	int id = 1; // Current individual to check
	double mindist;
	double newdist;

	while((ct < *N) && (id < *maxid))		// NOT <=*maxid because id starts at 0
	{ 
		mindist = 500;				// Smallest distance to any other point included in the sample
		for(int j = 0; j < id ; j++)		// Searches over all PREVIOUS points
		{
			if(use[new_order[j]])					// Calculates distance to points included in the sample already
			{
				newdist = rdist(x[new_order[j]], y[new_order[j]], x[new_order[id]], y[new_order[id]]);
				if( newdist < mindist ) // Keeps track of the smallest distance with any other point.
					mindist = newdist;
			}
		}

		if(mindist >= *buffer) //Point is only included if nearest neighbor is outside of buffer distance.
		{
			ct++;
			use[new_order[id]] = 1;
		}

		id++;
	}
	*N = ct; //Replaces # desired points with # actually fit.
}
 )