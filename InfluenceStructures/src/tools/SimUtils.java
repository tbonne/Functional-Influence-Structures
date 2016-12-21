package tools;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.collections.IteratorUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.geotools.referencing.GeodeticCalculator;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;

import LHP.Baboon;
import LHP.ModelSetup;
import LHP.Primate;
import LHP.Cell;
import LHP.Parameter;

import repast.simphony.space.gis.Geography;

public final strictfp class SimUtils {

	/* Method list:
	 * 		(1)	getCellAtThisCoord()
	 * 		(2) copyIterable()    				- converts Iterable to arrayList: primates
	 * 		(3) copyIterable()					- converts Iterable to arrayList: cells
	 * 		(4) getObjectsWithin()
	 */

	//(1)

	//(2)
	public static <T> ArrayList<Primate> copyIterable(Iterable<T> iter) {
		ArrayList<Primate> copy = new ArrayList<Primate>();
		while (iter.iterator().hasNext())
			copy.add((Primate) iter.iterator().next());
		return copy;
	}

	//(3)
	public static <T> ArrayList<Cell> copyIterable(Iterable<T> iter, boolean e) {
		ArrayList<Cell> copy = new ArrayList<Cell>();
		while (iter.iterator().hasNext())
			copy.add((Cell) iter.iterator().next());
		return copy;
	}

	//(4)
	public static  <T> List<T> getObjectsWithin(Coordinate coord,double radius, Class clazz){
		Envelope envelope = new Envelope();
		envelope.init(coord);
		envelope.expandBy(radius);
		Iterable<T> objectsInArea = ModelSetup.getGeog().getObjectsWithin(envelope,clazz);
		envelope.setToNull();

		//if(objectsInArea.iterator().hasNext()==false)System.out.println("warning no objects found of class: "+ clazz.toString());
		return IteratorUtils.toList(objectsInArea.iterator());
	}


	//find all food sites within a visual field
	/*public static  List<Cell> getFoodInFrontOfMe(Primate primate,double radius){

		ArrayList<Cell> inFront = new ArrayList<Cell>();
		double[] visAngles = {AngleUtils.getInitialAngle(primate.getFacing()),AngleUtils.getFinalAngle(primate.getFacing())};
		List<Cell> possiblePatches = SimUtils.getObjectsWithin(primate.getCoord(), Parameter.visual_range, Cell.class); 

		//radius; //distance
		for(Cell c : possiblePatches){

			double dist = c.getCoord().distance(primate.getCoord());

			//if within body radius
			if(dist<0.1){
				if(dist>0){
					if(isBlocked(primate,c.getCoord())==false){
						inFront.add(c);
					}
				}
			} else if(dist<radius){

				//less than left/right visual arc
				double angleToC = AngleUtils.getAngle(primate.getCoord(), c.getCoord());
				if(AngleUtils.withinAngles(primate.getFacing(), angleToC) ){
					if(isBlocked(primate,c.getCoord())==false){
						inFront.add(c);
					} 
				}
			} else if (Double.isNaN(dist)){

				if(isBlocked(primate,c.getCoord())==false){
					inFront.add(c);
				}
			}
		}

		if(inFront==null)System.out.println("no food in sight");

		return inFront;
	}*/

	public static List<Cell> getFoodAroundMe(Primate primate,double radius){
		ArrayList<Cell> nearby = new ArrayList<Cell>();
		List<Cell> possiblePatches = SimUtils.getObjectsWithin(primate.getCoord(), Parameter.visual_range, Cell.class); 

		//distance
		for(Cell c : possiblePatches){

			double dist = c.getCoord().distance(primate.getCoord());

			if(dist<radius){
				nearby.add(c);
			} 
		}
		return nearby;
	} 

	//find all primates within a visual field
	/*public static  List<Primate> getIndividualsInFrontOfMe(Primate primate,double radius){

		ArrayList<Primate> inFront = new ArrayList<Primate>();
		//GeodeticCalculator gc = ModelSetup.getGC();
		List<Primate> possiblePatches = SimUtils.getObjectsWithin(primate.getCoord(), Parameter.visual_range, Baboon.class);

		//radius; //distance
		for(Primate b : possiblePatches){
			if(b.getId()!=primate.getId()){
				//if individual is within visual distance
				if(b.getCoord().distance(primate.getCoord())<radius){

					//setup geo calculator
					//gc.setStartingGeographicPoint(primate.getCoord().x,primate.getCoord().y);
					//gc.setDestinationGeographicPoint(c.getCoord().x, c.getCoord().y);

					//less than left/right visual arc
					//if(Math.abs(primate.getFacing()-getAngle(primate.getCoord(),b.getCoord()))<Parameter.visualArc/2.0){ 
					if(AngleUtils.withinAngles(primate.getFacing(),AngleUtils.getAngle(primate.getCoord(),b.getCoord()))){
						inFront.add(b);
					}
				}
			}
		}
		return inFront;
	}


	private static boolean isBlocked(Primate primate,Coordinate c){

		boolean retval = false;

		//get direction to cell + distance
		double dist_to_food = primate.getCoord().distance(c);
		RealVector v_to_food = getVector(primate.getCoord(),c);
		v_to_food = v_to_food.unitVector();
		double dist_to_groupmate = Parameter.visual_range+1;
		RealVector v_to_groupmate = null;

		//get direction and distance to all other individuals: determine in one is in "front" of another
		for(Primate groupmate : ModelSetup.getAllPrimateAgentsInOrder()){
			if(primate.getId()!=groupmate.getId()){
				dist_to_groupmate = primate.getCoord().distance(groupmate.getCoord());
				v_to_groupmate = getVector(primate.getCoord(),groupmate.getCoord());
				v_to_groupmate = v_to_groupmate.unitVector();

				double vect_dist = v_to_groupmate.getDistance(v_to_food);
				if(vect_dist < Parameter.frontBlock_threshold && dist_to_food>dist_to_groupmate){
					retval=true;
				}
			}
		}

		return retval;

	}*/

	private static RealVector getVector(Coordinate source, Coordinate target) {
		//get coordinates from simple feature
		RealVector rv = new ArrayRealVector(new double[] { target.x - source.x, target.y - source.y }, false);
		return rv;
	}

	public static Coordinate generateCoordAround(double cenx,double ceny){
		double angle = Math.random()*Math.PI*2;
		double rDistance = Math.random()*Parameter.maxInitialGroupDistance;
		Coordinate testCoord = new Coordinate(cenx+Math.cos(angle)*rDistance,ceny+Math.sin(angle)*rDistance);

		while( (testCoord.x>0 && testCoord.y>0 && testCoord.x<Parameter.landscapeWidth && testCoord.y<Parameter.landscapeWidth) == false){
			angle = Math.random()*Math.PI*2;
			rDistance = Math.random()*Parameter.maxInitialGroupDistance;
			testCoord = new Coordinate(cenx+Math.cos(angle)*rDistance,ceny+Math.sin(angle)*rDistance);
		}

		return testCoord;

	}
	
	public static boolean occupied(Coordinate coord, int id){
		boolean overlap = false;
		for(Primate p:ModelSetup.getAllPrimateAgents()){
			if(p.getId()!=id){
				double dist = p.getCoord().distance(coord);
				if(dist<=Parameter.foodBuffer){
					overlap = true;
				}
			}
		}
		return overlap;
	}

}
