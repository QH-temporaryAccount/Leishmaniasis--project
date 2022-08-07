/**
* Name: Finalbatch
* Tags: Leishmania , Rodent , human, Transmission 
*/ model Finalbatch

global {
	list<float> TIME;
	list<int> INFECT;
	list<float> Temperatures;
	float current_temp;
	date starting_date <- date([2013, 1, 1]); // as not all months have the same number of days, it is better to use a real calendar. 
	float C_mosquito <- 0.02;
	float rodent_infection_prob <- 0.009;
 //Simulation starts on Januray 1st 2013.
 int
	current_month <- starting_date.month; // Dorra  
 // nb_infectes_step << roden count each.isInfected;
 // nb_sains_step << roden count not each.isInfected;
 reflex save_results
	when: MODE_BATCH {
		nb_infected_step << roden count each.isInfected;
		nb_notinfected_step << roden count not each.isInfected;
	}

	reflex manage_Temperatures {
		if (current_date.month != current_month) {
			int index <- months_between(starting_date, current_date);
			if (index < length(Temperatures)) {
				current_temp <- Temperatures[index];
				write "" + cycle + " date:" + current_date + " -> " + current_temp;
			} else {
				if (not MODE_BATCH) {
					do pause;
				}

			}

		}

		current_month <- current_date.month;
	} // Model Inputs files     
 file modis_data <- file("../includes/NDVIm1.tif");
	geometry shape <- envelope(modis_data);
	float MAX_FOOD <- 1.0;
	int dim_grille_width;
	int dim_grille_height;
	float SIZE_RODENT <- 10.0; //  Rodent size on the grid SIZE_RONGEUR
 // Parametres relative to the end of simulation  
 bool MODE_BATCH <- false;
	int SIMULATION_END <- 26280; //3*365*24;
 // paramettre relative to time and space scales  
 int CURRENT_HOUR update: (cycle / 60) mod 24;
	float step <- 1 #h; // Parametre relative to the disease   
 float Proba_infection <- 0.01;//1/100
	float alpha <- 0.5; //  positive value
 float CURING_PROBA<-0.01; //0.3/30
 
 
 
 
 
 
 
 
 
 
 
 float hh(float x){
 	float alpha1<-1.1;
 	return x^alpha1/(1+x^alpha1);
 }

 
 //parameters relative to human  unit:#day
 float CURING_PROBA_HUMAN<-1/(3*30);
 float infectionFactorHuman<-0.1;
 float miu<-0.01; //growth and death rate
 float incubation<-1/(5*7);
 float clinical_proportion<-2/3;
 
  float N<-SH+EH+IH+AH+RH;//human population
  float ratioMH<-sum(terrain_cell collect (each.IM)) / sum(terrain_cell collect (each.x)) update: sum(terrain_cell collect (each.IM)) / sum(terrain_cell collect (each.x));
  float IH<-0.0;//clinical compartment who get infected and show symptom
  float IH_stepPlusOne <- 0.0;
  float EH<-0.0;
  float EH_stepPlusOne<-0.0;
  float AH<-0.0;//Asymptomatic compartment
  float AH_stepPlusOne<-0.0;
  float RH<-0.0;
  float RH_stepPlusOne<-0.0;
  float SH<-1000.0;
  float SH_stepPlusOne<-0.0;

  
  reflex human_dynamic when: every(1 #day){
  float lambda<-infectionFactorHuman *hh(ratioMH); //force of infection for human
  
  	SH_stepPlusOne<-miu *N - lambda*SH-miu*SH;
 	SH<-SH_stepPlusOne +SH;
 	
 	EH_stepPlusOne<-infectionFactorHuman *hh(ratioMH)*SH-incubation*EH-miu*EH;
 	EH<-EH_stepPlusOne+EH;
 	
 	IH_stepPlusOne <- incubation * clinical_proportion*EH-CURING_PROBA_HUMAN*IH -miu*IH;
	IH <- IH_stepPlusOne+ IH;
	
	AH_stepPlusOne<-incubation * (1-clinical_proportion)*EH-CURING_PROBA_HUMAN*AH -miu*AH;
	AH <- AH_stepPlusOne+ AH;

	RH_stepPlusOne<-CURING_PROBA_HUMAN*(IH+AH)-miu*RH ;
	RH<-RH_stepPlusOne+ RH;
	
	//write "population: "+N;
	//write"SH: "+SH;
	//write"RH: "+ RH;
	
 }
 
 // parameters relative to mosquito
 float infectionFactor<-0.1;
 float mosquitoDeath<-0.002;//1/(3*7*24)
 
 // Parametres relative to rondent  
 int RODENT_NUM <- 50; // maximum number fo rodents 
 int
	nb_infected_init <- 1;
	float NEIGHB_INF_RADIUS <- 0.0;
	float FOOD_COMSUM <- 0.001;
	float MIN_FOOD_DETECTION <- 0.7; // the minimum amount of resource that a rodent can see, 
 // below this value the rodent cannot 

	// see the resource and therefore cannot eat it
 float LOST_MEM_RATE <- 0.2; // Loose memory rate per day
 // parametre relative to temperature and precipitation

	// float Temperature <- 7 * sin(360 * (time / (365 * 24 * 3600))) + 20 update: 7 * sin(360 * (time / (365 * 24 * 3600))) + 20;

	/*float Precipitation <- 27.5*sin(2°pi*(time/(365*12)))+5 update:27.5*sin(2°pi*(t
e/(365*12)))+5;  */ 
	int cpt_ndvi <- 0; /***OUTPUT DATA AND PARAMETERS * */
	// output relative to reodents
 int Prev_rodent_num_infected <- nb_infected_init update: Rodent_num_infected;
	int Rodent_num_infected <- nb_infected_init;
	float Infected_rate update: empty(roden) ? 0 : Rodent_num_infected / length(roden);
	int Rodent_num_not_infected <- (RODENT_NUM - Rodent_num_infected) update: Rodent_num_infected;
	int Rodent_num_newinfected <- 0;
	float R_0 <- 0.0;
	float Mean_food_consom; // Mean consuming food by rodent
 float Mean_rand_mv_num; //  Mean number of random movement  number of #rondent 0
 float Mean_pref_mv_num;
	// Mean number of preferential  movement  number of #rondent 0
 list<int> nb_infectes_step;
	list<int> nb_sains_step;
	list<float> Mean_rand_mv_num_step;
	list<float> Mean_pref_mv_num_step; //output relative to the environment 
 float Mean_env_food; //Mean food in the environement
 //reflex infected_number  when: every(1#h)  {

	// save the values of the variables name, speed and size to the csv file; the rewrite facet is set to true to continue to write in the same file

	//save [Rodent_num_infected] to: "../results/seed35.csv" type:"csv" rewrite: false;
 //Pause the model as the data are saved
 //}
 reflex pausing when: not MODE_BATCH and
	(cycle = SIMULATION_END or Rodent_num_infected = 0) {
		do pause;
	}

	reflex save_data when: false and every(1 #h) {
		Rodent_num_infected <- roden count (each.isInfected);
		Rodent_num_not_infected <- RODENT_NUM - Rodent_num_infected;
		Rodent_num_newinfected <- Rodent_num_infected - Prev_rodent_num_infected;
		R_0 <- Rodent_num_not_infected != 0 ? (Rodent_num_infected - Prev_rodent_num_infected) / (Rodent_num_not_infected) : 0;
		Mean_food_consom <- mean(roden collect each.foodEated);
		Mean_env_food <- mean(terrain_cell collect each.food);
		Mean_rand_mv_num <- mean(roden collect (each.randMv / (each.randMv + each.prefMv)));
		Mean_pref_mv_num <- mean(roden collect (each.prefMv / (each.randMv + each.prefMv))); // save [Rodent_num_infected] to: "../results/seed6.csv" type:"csv" rewrite: false;

		//Pause the model as the data are saved
 }

	reflex pausing when: not MODE_BATCH and (cycle = SIMULATION_END) {
		do pause;
	} // Saving  data and  stop simulation 
 reflex save_results when: MODE_BATCH {
		nb_infectes_step << roden count each.isInfected;
		nb_sains_step << roden count not each.isInfected;
		Mean_rand_mv_num_step << mean(roden collect (each.randMv / (each.randMv + each.prefMv)));
		Mean_pref_mv_num_step << mean(roden collect (each.prefMv / (each.randMv + each.prefMv)));
	}

	reflex update_cell when: every(1 #month) {
		int id <- 2 + (cpt_ndvi mod 15);
		string modis_file_path <- ("../includes/NDVIm" + id + ".tif");
		write "" + cycle + " -> " + file_exists(modis_file_path);
		if (file_exists(modis_file_path)) {
			grid_file data <- (grid_file(modis_file_path));
			ask terrain_cell {
				food <- float(data[grid_y * dim_grille_width + grid_x] get "grid_value");
			}

		}

		cpt_ndvi <- cpt_ndvi + 1;
	}

	reflex update_memory_grid { //    
 roden rod <- roden(0);
		matrix<float> data1 <- rod.infoMatrix;
		matrix<float> data2 <- rod.foodMatrix;
		ask memory_grid {
			memory <- float(100 * data1[grid_x, grid_y] * data2[grid_x, grid_y]);
		}

	} // Initialization         
 init {
		file my_csv_file <- csv_file("../includes/Temp2.csv", ",");
		loop l over: my_csv_file.contents {
			Temperatures << float(l);
		}

		current_temp <- Temperatures[0];
		dim_grille_height <- matrix(terrain_cell).rows;
		dim_grille_width <- matrix(terrain_cell).columns;
		ask terrain_cell {
			food <- 1 - (((color as list) at 0) / 255);
		} // Agent creation       
 create roden number: RODENT_NUM;
		ask nb_infected_init among roden {
			isInfected <- true;
		}

	}

}


/*Define the terrain: a grid of square cells with dim_grid sides and each cell has 8 neighboring cells. 
 Then  define a plant density in each cell and therefore the color of the cell will depend on the plant density*/
grid terrain_cell file: modis_data neighbors: 8 { //file: MODIS_DATA{
 // Grid representing the environment with food proportion to ndvi of satelite fig

//  list<float> mouvementProb; // list of the probability of moving to another cell
 float foodProd <- rnd(10) / 10; /*Precipitation/50 ;*/ 
 float food <- grid_value min: 0.0 max:
	MAX_FOOD;
	rgb color <- rgb(int(255 * (1 - food)), 255, int(255 * (1 - food))) update: rgb(int(255 * (1 - food)), 255, int(255 * (1 - food)));
	float h <- 100.0;
	float t; // time 
 float x <- 4.0; // mosquito numbers 
 float Temp <- current_temp / 30.36; // temp scaling , maxmium temp is 30.36
 float IM <- 0.01;
	float SM; 

//carrying capacity-k , updated by NDVI- food.cell-MAX_FOOD, logistic
float max_K<-500.0;
	float k update: max_K*food/(1+food); // allee function
	
 // gauss for pp , updated by NDVI and temperture, 2 dimensions
 float ppFunction (float x1, float x2, float
	x3, float x4) {
		return normal_density(x1, x2, x3) + x4 / (1 + x4);
	}
	float tempOptimal <- 27.0;
	float tempMax <- 30.36;
	float tempMin <- 10.55;
	float tempSig <- sqrt((tempOptimal - tempMax) ^ 2 + (tempOptimal - tempMin) ^ 2) / 2;
	float pp update: ppFunction(current_temp, tempOptimal, tempSig, food); 
 // just to test
 list<roden> RodentHereList update: roden inside (self);
	//list<roden> RodentHereList update: roden at_distance (1) ;
 int infectedHere <- RodentHereList count (each.isInfected = true) update: RodentHereList count
	(each.isInfected = true); //  = infected rodent in the patch 
 int notinfectedHere <- RodentHereList count (each.isInfected = false) update: RodentHereList count
	(each.isInfected = false); //  = infected rodent in the patch 
 float ratio update: empty(RodentHereList) ? 0 : infectedHere / length(RodentHereList);

	/*  equation eqX {
    diff(self.x,t) = 0.02 * pp * self.x * (1 - self.x / (400)) + self.x;
   // diff(x,t) =  1/exp(x) ;
    diff(IM,t) =  (x- IM) * ratio;
  //  diff(IM,t) =  infectionFactor * ratio * (x - IM);
    } // to be checked */
	float x_stepPlusOne <- 0.0;
	float IM_stepPlusOne <- 0.0; //Mosquitoes Dynamic
 reflex Verhulst {
		x_stepPlusOne <- C_mosquito * pp * x * (1 - x / (k)) + x;
		//pp : growth rate dep on temp and ndvi, (27temp, max indvi ) max number of mosquito in the scale carrying capacity in the cell : 400 function of the ndvi  , check pp add a scale factor // zoom factor

		x <- x_stepPlusOne;
	} //Infected Mosquitoes

	reflex MosquitoEpidemie {
		IM_stepPlusOne <- infectionFactor * ratio * (x - IM) - mosquitoDeath*IM + IM;
		IM <- IM_stepPlusOne; //      write " IM " + self.IM +  " ratio" + ratio + "   x   " + self.x + "  old IM"  + self.IM  ;
 }
  
	/* 
reflex solving  when:  every(1 # day){
	    write "pp =  " + pp ;
		solve eqX method: "rk4" step_size: h;
}*/ }

grid memory_grid file: modis_data neighbors: 8 { // Grid representing the memory of rodent #0   
 float memory;
	rgb color <- rgb(int(255 * (1 - memory)), 255, int(255 * (1 - memory))) update: rgb(int(255 * (1 - memory)), 255, int(255 * (1 - memory)));
} /*Rodents move randomly and can consume resources */ 

species roden skills: [moving] { // Main agent of the model
 list<roden> neighbInfectedRodent update: roden at_distance
	(NEIGHB_INF_RADIUS) where (each.isInfected = false); // non infected rodent in the neighb 
 /*terrain_cell myCell <- one_of(terrain_cell);*/ float curingProb;
	float infectionProb;
	float foodEated;
	bool isInfected;
	float minFoodDetection <- MIN_FOOD_DETECTION;
	int prefMv; // counter to know the number of prefered mv 
 int randMv; // counter to know the number of random  mv
 terrain_cell myCell <- one_of(terrain_cell[1, 1].neighbors);
	/*Initialization of rodent location */ matrix<float> infoMatrix; /* Info matrix  */ matrix<float> foodMatrix; /* Estimation matrix */
	/***********Initialisation des agents *************** */ init {
		prefMv <- 1;
		randMv <- 1;
		infectionProb <- Proba_infection; // 
 curingProb <- CURING_PROBA;
		foodEated <- 0.0;
		infoMatrix <- 0.0 as_matrix ({dim_grille_width, dim_grille_height});
		foodMatrix <- 0.0 as_matrix ({dim_grille_width, dim_grille_height}); //location<-myCell.location;
 } // Agetns shape
 aspect base {
		draw circle(SIZE_RODENT) color: isInfected ? #red : #blue;
	}

	aspect info {
		draw circle(SIZE_RODENT) color: color;
		draw string([myCell.grid_y + 1, myCell.grid_x + 1]) size: 3 color: #black;
	} /************Prcedure des agents WAJDY *************************** */ //    action calcule_info (terrain_cell cell){

	//        float sumc <- 0.0; /*somme des cellules connues */
 //        float sumr <- 0.0; /*somme des cellules ayant des ressources non nulles */

	//        int count <- 0; /*nombre de cellules voisines */
 //        float a <- 0.0;
 //        loop i over: cell.neighbors{

	//            sumc <- M[i.grid_x, i.grid_y] != -1.0 ? sumc + 1 : sumc; /*Calcule de la somme des cellules connues */

	//            sumr <- M[i.grid_x, i.grid_y] = 1.0 ? sumr + 1 : sumr; /* Calcule de la somme des cellules ayant des ressources */

	//            count <- count + 1; /*Nombre de cellules voisines */
 //        }
 //
 //        infoMatrix[cell.grid_x, cell.grid_y] <- sumc / (count);

	//        a <- infoMatrix[cell.grid_x, cell.grid_y];
 //        foodMatrix[cell.grid_x, cell.grid_y] <- (sumr / sumc); //(sumr / count) * (1 / a);
 //    }
 //

	//    action collecte_info (terrain_cell cell){

	//    /*Pour toute cellule voisine on va calculer les matrices d'information, de certitude et de ressource sachant ses voisines */

	//        M[cell.grid_x, cell.grid_y] <- cell.food > 0.3 ? 1 : 0;
 //        infoMatrix[cell.grid_x, cell.grid_y] <- M[cell.grid_x, cell.grid_y];

	//        loop i over: cell.neighbors
 //        {
 //            do calcule_info(i);
 //        }
 //
 //    }
 //    
 action matrix_update (terrain_cell cell) {
	/*For any neighboring cell we will compute the information, info/certainty and resource matrices giving its neighbors */ // Info Collection 

		foodMatrix[cell.grid_x, cell.grid_y] <- cell.food > MIN_FOOD_DETECTION ? cell.food : 0;
		infoMatrix[cell.grid_x, cell.grid_y] <- 1.0; //M[cell.grid_x, cell.grid_y];
 matrix<float> Minfo <- copy(infoMatrix);
		matrix<float> Mres <- copy(foodMatrix); // Filling Matrices 
 loop i over: cell.neighbors {
			float sumI <- 0.0; /*sum of known cells */ float sumR <- 0.0; /*sum of cells with non-zero resources */ loop k over: i.neighbors {
				sumI <- infoMatrix[k.grid_x, k.grid_y] + sumI;
				sumR <- foodMatrix[k.grid_x, k.grid_y] + sumR;
			}

			int count <- length(i.neighbors);
			Minfo[i.grid_x, i.grid_y] <- max([infoMatrix[i.grid_x, i.grid_y], sumI / count]); //sumI/count
 Mres[i.grid_x, i.grid_y] <- infoMatrix[i.grid_x, i.grid_y] != 1 ?
			sumR / count : foodMatrix[i.grid_x, i.grid_y];
		}

		foodMatrix <- copy(Mres);
		infoMatrix <- copy(Minfo);
	}

	action cell_to_go (terrain_cell cell) { // Choice of the destination cell with respect to the foodMatrix and infoMatrix

	// The chosen cell is the one that maximizes the forward product of the two matrices in the rodent's neighborhood

	// If the max is lower than a certain value the movement is random
 float minValue <- 10 ^ (-4);
		float a <- infoMatrix[cell.grid_x, cell.grid_y];
		float b <- foodMatrix[cell.grid_x, cell.grid_y];
		float max <- a * b;
		terrain_cell destination_cell <- cell;
		loop i over: cell.neighbors {
			a <- infoMatrix[i.grid_x, i.grid_y];
			b <- foodMatrix[i.grid_x, i.grid_y];
			if a * b > max {
				max <- a * b;
				destination_cell <- i;
			}

		}

		if max < minValue {
			destination_cell <- one_of(cell.neighbors);
			randMv <- randMv + 1;
		} else {
			prefMv <- prefMv + 1;
		}

		return destination_cell;
	}

	reflex lost_memory { //} when:  every(1 # day){
 infoMatrix <- infoMatrix * (1 - LOST_MEM_RATE);
	}

	reflex move_eat { // Moving and matrix update
 terrain_cell dest;
		if (myCell.food < MIN_FOOD_DETECTION) {
			dest <- cell_to_go(terrain_cell(self.location)); //myCell);
 do goto target: dest.location speed: 1.0;
			myCell <- dest;
			do matrix_update(myCell);
		} else {
			myCell.food <- myCell.food - min([myCell.food, FOOD_COMSUM]);
			foodEated <- min([myCell.food, FOOD_COMSUM]);
		}

	} //    }
 /********************Epidemiology*************************************************************
 */ float puissance (float a, float b, float x, float y) {
	// a and b are both positive
 // float b=normal_density(2,5,1)
 a <- 10 ^ (-5);
		b <- 10 ^ (-5);
		return x ^ a * y ^ b;
	}

	float som (float a, float x, float y) { // float b=normal_density(2,5,1)
 return a * x + (1 - a) * y;
	}
	/*   	float gauss (float x, float y){
   
         float alp<-1.5; // parameter to fit the effect of temperature 
         float tempMean<-19.95;
         float tempMax<-30.36;
         float tempMin<-10.55;
         float sig<-sqrt((tempMean-tempMax)^2+(tempMean-tempMin)^2)/2;
         return (alp*y/(sig*sqrt(2*#pi))*exp(-0.5*((x-tempMean)/sig))^2) ;
         }
*/
	float f (float x) {
		return x / (1 + x);
	}

	reflex infect when: (myCell.IM != 0) {
		float Temp <- current_temp / 30.36; //     Tempature scaling  
 list<roden> rodentHere <- roden inside (myCell);
		float ratioMR <- 1.0;
		if (!empty(rodentHere)) {
			ratioMR <- f(myCell.IM / length(rodentHere)); // add a function of the ratio; to be defined x/1+x 
 }

		if flip(ratioMR * rodent_infection_prob) { // infection prob  global 
 isInfected <- true;
		}

	}

	reflex guerison when: isInfected and every(1 #week) {
		if flip(CURING_PROBA) {
			isInfected <- false;
		}

	} /*************************************************************** */ }

experiment simulation type: gui {

	float seedValue <- 1.0;
	float seed <- seedValue; // forces the value of the seed.
	
	parameter "infectionFactor" var:infectionFactor <-0.55;//fit mosquito infectionFactor to mosquito ratio data

//parameters for sensitivity analysis: RODENT_NUM and max_K
  parameter "RODENT_NUM" var:RODENT_NUM<-50;
  
 output {
 
// monitor "distance_variation" value: new_distance / 30;
//monitor "nb_people_newinfected" value: Rodent_num_newinfected;
//monitor "R_0" value: R_0; 
//monitor "occupied cells" value: terrain_cell count (each.is_occupied = true);
// monitor "Total distance" value: total_distance;

		/*monitor"Temperature" value: Temperature;*/
		/*display map  {
            species roden aspect:base;
            grid terrain_cell lines: #black;            
        }
        * 
        */ 
        
        
        
        
        display main_display {
			grid terrain_cell lines: #black; // size: { 1, 0.5 } position: { 0, 0 };
 species roden aspect: base;
			//            grid toto lines: # black size: { 1, 0.5 } position: { 0, 0.5 };
 //           species roden aspect: base;
 }

		/*display memoire refresh: every(5 #cycle) {
			grid memory_grid lines: #black; // size: { 1, 0.5 } position: { 0, 0 };
             
             grid toto lines: # black size: { 1, 0.5 } position: { 0, 0.5 };

			           species roden aspect: base;
 //        
 }*/

		/*display Temperature_information refresh: every(5 #cycle) {
			chart "Temperature" type: series size: {1, 0.5} position: {0, 0} {
				data "Temperature" value: current_temp color: #red;
			}

			chart "Food evolution" type: series size: {1, 0.5} position: {0, 0.50} {
				data "food" value: terrain_cell count (each.food > MIN_FOOD_DETECTION) color: #blue;
			}

			chart "Food evolution" type: series size: {1, 0.5} position: {0, 0.50} {
				data "Mean Consom food" value: (100 * Mean_food_consom) color: #red;
				data "Mean food" value: Mean_env_food color: #green;
			}

		}*/

		/*display mouvement_information refresh: every(5 #cycle) {
			chart "mouvement " type: series {
				data "random mv" value: Mean_rand_mv_num color: #green;
				data "priv mv" value: Mean_pref_mv_num color: #red;
			}

		}*/

		display Infection_information refresh: every(5 #cycle) {
			chart "Disease spreading" type: series size: {1, 0.5} position: {0, 0} {
				//data "susceptible" value: roden count (each.isInfected = false) color: #green;
				data "infected" value: roden count (each.isInfected = true) color: #red;
				data "infected Human" value: IH color:#blue;
				data "infected Mosquito" value: sum(terrain_cell collect (each.IM)) color: #black;
			}

		}

		display human refresh: every(1 #day) {
			chart "Disease spreading" type: series size: {1, 0.5} position: {0, 0} {
				//data "infected" value: mean(terrain_cell collect (each.IM / each.x)) color: #red;
				data "susceptible" value:SH color:#black;
				data "exposure" value:EH color:#green;
				data "asymptomatic" value:AH color:#yellow;
				data "infected" value: IH color:#red;
				data "recovery" value: RH color:#blue;
				
			}

		}

		display ratios refresh: every(5 #cycle) {
			chart "Mosquito Disease spreading" type: series size: {1, 0.5} position: {0, 0} {
				data "ratio Mosquito" value: sum(terrain_cell collect (each.IM)) / sum(terrain_cell collect (each.x)) color: #red;
				data "ratio rodent " value: roden count (each.isInfected = true) / length(roden) color: #orange;
				data "human ratio" value: IH/N color:#green;
				}
			}
			
			display mosquitoDynamic refresh: every(5 #cycle){
				chart "Mosquito Disease number" type: series size: {1, 0.5} position: {0, 0} {
			      data "total mosquito" value: sum(terrain_cell collect (each.x)) color:#blue;
			      data "infected mosquito" value: sum(terrain_cell collect (each.IM)) color:#red;
				}
			}
}



	
	reflex save_infect when:  every(5 #cycle) {
		  float mosquitoInfectionRatio <- sum(terrain_cell collect (each.IM)) / sum(terrain_cell collect (each.x));
		  float rodentInfectionRatio <- roden count (each.isInfected = true) / length(roden);
		  float humanInfectionRatio <- IH/N;
		  float humanSusceptibleNumber<-SH;
		  float humanExposureNumber<-EH;
		  float humanInfectNumber<-IH;
		  float humanRecoverNumber<-RH;
		  float humanAsymptomaticNumber<-AH;
		  float mosquitoDensity<-sum(terrain_cell collect (each.x));
		  float mosquitoInfectNumber<-sum(terrain_cell collect (each.IM));
		  
			// save the values of the ratios, human dynamic, mosquito dynamic to the csv file; the rewrite facet is set to false to continue to write in the same file
		  save [mosquitoInfectionRatio, rodentInfectionRatio, humanInfectionRatio] to: "../includes/ratio.csv" type:"csv" rewrite: false;
		  save [humanSusceptibleNumber, humanExposureNumber, humanAsymptomaticNumber, humanInfectNumber, humanRecoverNumber] to: "../includes/human_dynamic.csv" type:"csv" rewrite: false;
		  save [mosquitoDensity, mosquitoInfectNumber] to: "../includes/mosquito_dynamic.csv" type:"csv" rewrite: false;
		  }
		  }


