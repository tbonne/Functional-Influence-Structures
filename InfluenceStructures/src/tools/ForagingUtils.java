package tools;

import LHP.Cell;
import LHP.ModelSetup;
import LHP.Parameter;
import LHP.Primate;
import LHP.Baboon;
import repast.simphony.random.RandomHelper;

import com.vividsolutions.jts.geom.Coordinate;


public final strictfp class ForagingUtils {

	/* Method list:
	 * 		(1)	choosefeedingSite()							- compares visual and remembered sites, chooses best one
	 * 		(2) findBestUnknownFoodSource()					- finds best visual site
	 * 		(3) findBestKnownFoodSourceTopo()				- finds best remembered site
	 * 		(4)	chooseFoodSourceTowardsMyMigrationGoal()	- finds best site on way to my chosen site
	 */

/*
	//(1)
	//finds the best resource
	public static Coordinate findBestUnknownFoodSource(Primate primate){//find the best unknown resource; all resources within search radius

		//if(primate.getId()==13){
		//	System.out.println("choosing unknown site");
		//}
		
		
		//System.out.println("choosing unknown site");
		double dist_I_F = 0;
		double resourceWeight = 0;
		double minResourceWeight = 999999;
		Coordinate bestCloseSite = null;

		//Calculate each resources food weight and choose the lowest one (weight = distance / food)
		for (Cell r:primate.getVisualPatches()){

			//if(closest(primate, r)==true){
			dist_I_F = primate.getCoord().distance(r.getCoord());
			//dist_I_F = primate.getBaboonFollower().getCoord().distance(r.getCoord());
			//dist_I_F = Math.max(dist_I_F, Parameter.bodyRadius);

			//within 4m of another?
			boolean close = true;
			for(Primate p : ModelSetup.primatesAll){
				if(p.getId()!=primate.getId()){
					double di = p.getCoord().distance(r.getCoord());
					if (di<Parameter.bodyRadius && di<dist_I_F)close=false;
				}
			}

			if(close == true){
				resourceWeight=dist_I_F/r.getResourceLevel();
				if(resourceWeight<minResourceWeight){
					minResourceWeight = resourceWeight;
					bestCloseSite = r.getCoord();
				}
			}
		}

		if(bestCloseSite==null){
			//System.out.println("no food site within view " +primate.getId());
			//Coordinate rand = new Coordinate(primate.getCoord().x+RandomHelper.nextDoubleFromTo(-Parameter.visual_range/5, Parameter.visual_range/5),primate.getCoord().y+RandomHelper.nextDoubleFromTo(-Parameter.visual_range/5, Parameter.visual_range/5));
			primate.setFacing(primate.getFacing()+(Math.random()-0.5)*40);
			Coordinate rand = new Coordinate(primate.getCoord().x+Math.cos(Math.toRadians(primate.getFacing()))*Parameter.maxDistancePerStep,primate.getCoord().y+Math.sin(Math.toRadians(primate.getFacing()))*Parameter.maxDistancePerStep);
			bestCloseSite = rand;
			//primate.setFoodTarget(rand);
			//look for a way around
			/*int i=0;
			if(RandomHelper.nextDouble()>0.5){
				i=-1;
			} else {
				i=1;
			}
			primate.setFacing(primate.getFacing()+(10*i));///
		}

		return bestCloseSite;
	}
*/
	/*public static Coordinate findBestUnknownFoodSource(Primate primate, Coordinate rand){//find the best unknown resource; all resources within search radius

		//System.out.println("choosing unknown site");
		double dist_I_F = 0;
		double resourceWeight = 0;
		double minResourceWeight = 999999;
		Coordinate bestCloseSite = null;

		//Calculate each resources food weight and choose the lowest one (weight = distance / food)
		for (Cell r:primate.getVisualPatches()){

			//if(closest(primate, r)==true){
			dist_I_F = primate.getCoord().distance(r.getCoord());

			//within 4m of another?
			boolean close = true;
			for(Primate p : ModelSetup.primatesAll){
				double di = p.getCoord().distance(r.getCoord());
				if (di<Parameter.bodyRadius)close=false;
			}

			if(close == true){
				resourceWeight=dist_I_F/r.getResourceLevel();
				if(resourceWeight<minResourceWeight){
					minResourceWeight = resourceWeight;
					bestCloseSite = r.getCoord();
				}
			}
		}

		if(bestCloseSite==null){
			//System.out.println("no food site within view " +primate.getId());
			bestCloseSite = rand;
			//primate.setFoodTarget(rand);
		} else {
			primate.setFoodTarget(null);
		}

		return bestCloseSite;
	}

	 */
}
