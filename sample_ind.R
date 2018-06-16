cppFunction(
'LogicalVector new_sample(NumericMatrix xy, int N, double buffer, LogicalVector use, IntegerVector new_order)
{
  int maxid = use.size();
	double newdist;
  bool * trydist = new bool[maxid];
  buffer = buffer * 1000; // km to m
  
  int ct = 0;
  for(int i = 0; i < maxid; i++){
    ct = ct + use[i];
    trydist[i] = 1;   // initializes trydist
  }

int id;  
  if( ct == 0 ){     
	    use[new_order[0]] = 1; // Automatically include the first point in the sample = 1st wolverine
      ct++;                  // Number of individuals included so far
      
      id = 1;  // If starting from scratch, add first ind and start searching at second random individual
      for(int i = 0; i < maxid; i++)
      {
        newdist = sqrt(pow(xy(new_order[id],1) - xy(i,1), 2)+ pow(xy(new_order[id],0) - xy(i,0), 2));
        if(newdist <= buffer) 
          trydist[i] = 0;
      } 
      
	} else {   
    for(int j = 0; j < maxid; j++){      
      if(use[j] == 1){
        for(int i = 0; i < maxid; i++)
        {
          newdist = sqrt(pow(xy(j,1) - xy(i,1), 2)+ pow(xy(j,0) - xy(i,0), 2));
          if(newdist <= buffer) 
            trydist[i] = 0;
        }
      }
    }
  id = 0;    // If starting from a population, start searching at first random
  }    
 
 	while((ct < N) && (id < maxid))		// NOT <=*maxid because id starts at 0
	{ 
		if(trydist[new_order[id]] == 1){
      use[new_order[id]] = 1;
      ct++;
     for(int i=0; i < maxid; i++)
     {
      newdist = sqrt(pow(xy(new_order[id],1) - xy(i,1), 2)+ pow(xy(new_order[id],0) - xy(i,0), 2));
      if(newdist <= buffer) 
        trydist[i] = 0;    
     }
    } 
		id++;
	}
  delete [] trydist;
  return use;
 }'
 )