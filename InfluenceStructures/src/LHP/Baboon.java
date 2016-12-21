package LHP;

import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import repast.simphony.random.RandomHelper;
import tools.SimUtils;
import tools.MoveUtils;
import jsc.distributions.Lognormal;
import cern.jet.random.*;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;

public class Baboon extends Primate{

	/****************************
	 * 							*
	 * Building a red colobus 	*
	 * 							*
	 * *************************/

	//initialize a primate agent
	public Baboon(int ID, Coordinate c, int groupSize, boolean isMale, int group){
		this.id = ID;
		this.coordinate = c;
		this.myGroup=group;
		if (isMale == true){
			this.sex =1; 
		} else {
			this.sex = 0;
		}
		//groupMates = new ArrayList<Primate>();
		destination = null;
		RealVector initialFacing = new ArrayRealVector(2,0);
		initialFacing.addToEntry(0,RandomHelper.nextDoubleFromTo(-1, 1));
		initialFacing.addToEntry(1,RandomHelper.nextDoubleFromTo(-1, 1));
		//initialFacing.addToEntry(0,0*RandomHelper.nextDouble()*0.001);
		//initialFacing.addToEntry(1,-1*RandomHelper.nextDouble()*0.001);
		this.setFacing(initialFacing);
		blocked=false;
		foodTarget=null;
		followMate=null;

		//vm = Parameter.vm;
		myVector = new ArrayRealVector(2,0.1);
	}


	/************************************
	 * 									*
	 * Stimuli (internal + external) 	*
	 * 									*
	 * *********************************/

	public void getInputs(){
		//update where the agent is on the landscapes
		this.coordinate = ModelSetup.getAgentGeometry(this).getCoordinates()[0];
		
		//update which patch i'm in
		this.myPatch = this.getCurrentFoodPatch();
		//i'm here (updating the code to have network structure implicit in the design of influence... i.e., adding in empirically extracted values)
		//update which patches i can see
		if(this.myPatch!=null){
			this.visualPatches = this.getVisibleFoodPatches(Parameter.visual_range,myPatch.visibleSites);
		} else {
			this.visualPatches = this.getVisibleFoodPatches(Parameter.visual_range);
		}


		//update my social partner of interest (attraction)
		if(followMate==null){
			Collections.shuffle(this.primateList);
			for(Primate pp : this.primateList){
				if(pp.getMyGroup()==this.getMyGroup() && pp.getId()!=this.getId()){
					followMate=pp;
					break;
				}
			}
		} else if (Math.random()<Parameter.followMateProb){
			Collections.shuffle(this.primateList);
			for(Primate pp : this.primateList){
				if(pp.getMyGroup()==this.getMyGroup() && pp.getId()!=this.getId()){
					followMate=pp;
					break;
				}
			}
		}
	
	
		//move based on attraction to food/social partners and repulsion from strangers
		attractionRepulsionRand();
		
	}

	/****************************
	 * 							*
	 * Behavioural response 	*
	 * 							*
	 * *************************/

	public void behaviouralResponse(){
		if(myPatch!=null){
			if(myPatch.getResourceLevel()>Parameter.biteSize){
				myPatch.eatMe(Parameter.biteSize);
			} else {
				move(myVector,true);
			}
		} else {
			move(myVector,true);
		}
	}

	/****************************
	 * 							*
	 * 		Energy updates		*
	 * 							*
	 * *************************/

	public void energyUpdate(){
		//nothing right now: agents are assumed to be always hungry, and will eat when in a patch
	}

	/****************************
	 * 							*
	 * 	   Behavioural model	*
	 * 							*
	 * *************************/


	private void attractionRepulsionRand(){
		
		//// calculate adjustment for previous bearing
		myVector.unitize();
		myVector = myVector.mapMultiply(Parameter.bearingWeight);

		//// calculate adjustment for food environment
		//calculate the weights of each patch
		ArrayList<Double> weights = new ArrayList<Double>();
		ArrayList<RealVector> directions = new ArrayList<RealVector>();
		double sum = 0;

		for(Cell c : this.visualPatches){

			//calculate weight
			double distance = Math.max(Parameter.foodBuffer,  c.getCoord().distance(this.getCoord()));
			double weight = (c.getResourceLevel()) / distance;  //no effect of     *  c.getFamiliarity(this.getId() 
			weights.add(weight);
			sum = sum + weight;

			//calculate the direction
			RealVector patchVector = new ArrayRealVector(2);
			double distX = c.getCoord().x-this.getCoord().x;
			double distY = c.getCoord().y-this.getCoord().y;
			if(distX!=0 && distY!=0){
				patchVector.setEntry(0,distX);
				patchVector.setEntry(1,distY);
				patchVector.unitize();
			}
			directions.add(patchVector);
		}
		//standardize the patch weights to sum to one
		for(Double d : weights){
			d=d/sum;
		}
		//calculate the avg direction (weights*dir to each patch)
		RealVector avgFoodVector = new ArrayRealVector(2);
		boolean foodFound = false;
		for(int i = 0 ; i< weights.size();i++){
			avgFoodVector = avgFoodVector.add(directions.get(i).mapMultiply(weights.get(i)));
			foodFound = true;
		}
		//add to myVector
		if(foodFound==true){
			avgFoodVector.unitize();
			avgFoodVector = avgFoodVector.mapMultiply(Parameter.foodWeight);
			myVector = myVector.add(avgFoodVector);
		}
		
		////calculate adjustment for attraction
		double distToPartner = followMate.getCoord().distance(this.getCoord());
		double magnitudeAtt = Math.max(Parameter.attractionWeight*(1-(Parameter.attractionDistMax/distToPartner)),0);
		double angleToP = MoveUtils.getAngle(this.getCoord(), followMate.getCoord(),true);
		myVector.addToEntry(0, magnitudeAtt*Math.cos(angleToP));
		myVector.addToEntry(1, magnitudeAtt*Math.sin(angleToP));
		
		////calculate adjustment for repulsion
		double magnitudeRep = 0;
		for(Primate stranger:getStrangers()){
			double distToStr = stranger.getCoord().distance(this.getCoord());
			magnitudeRep = Math.max(0, Parameter.repulsionWeight*((1-(distToStr/Parameter.repulsionDistMax))));
			double angleToS = MoveUtils.getAngle(this.getCoord(), stranger.getCoord(),true);
			myVector.addToEntry(0, -magnitudeRep*Math.cos(angleToS));
			myVector.addToEntry(1, -magnitudeRep*Math.sin(angleToS));
		}
		//// Choose final vector based on uncertainty around alternative influencing factors
		double u = Math.atan2(myVector.getEntry(1), myVector.getEntry(0));
		double length = Math.pow(Math.pow(myVector.getEntry(0),2)+Math.pow(myVector.getEntry(1), 2),0.5);
		double maxLength = Parameter.bearingWeight + Parameter.foodWeight + magnitudeAtt + magnitudeRep;
		double k = Math.max(-2*Math.log(length/maxLength),0.000001);
		 if(k<=0){
			 System.out.println("zeero u");
			 k=0.0;
		 }
		if(u==0)System.out.println("zeero u");
		if(u!=0){
			u = u + VonMises.staticNextDouble(1/k); 
			double deltaX = Math.cos(u);
			double deltaY = Math.sin(u); 
			RealVector finalVector = new ArrayRealVector(2);
			finalVector.setEntry(0, deltaX);
			finalVector.setEntry(1, deltaY);
			finalVector.unitize();
			myVector = finalVector;
		}
	}


	/****************************
	 * 							*
	 * 	   Methods				*
	 * 							*
	 * *************************/


	private ArrayList<Primate> getStrangers(){
		ArrayList<Primate> strangers = new ArrayList<Primate>();
		for(Baboon bab : this.primateList){
			if(this.coordinate.distance(bab.coordinate)<Parameter.repulsionDistMax && bab.myGroup!=this.myGroup){
				strangers.add(bab);
			}
		}
		return strangers;
	}

	private ArrayList<Cell> getVisibleFoodPatches(double f){

		Iterable<Cell> objectsInArea = null;
		Envelope envelope = new Envelope();
		envelope.init(this.getCoord());
		envelope.expandBy(f);
		objectsInArea = ModelSetup.getGeog().getObjectsWithin(envelope,Cell.class);
		envelope.setToNull();

		ArrayList<Cell> obj = new ArrayList<Cell>();
		while(objectsInArea.iterator().hasNext()){
			Cell neigh = objectsInArea.iterator().next();
			if(neigh.getCoord().distance(this.getCoord())<f){
				boolean occupied = SimUtils.occupied(neigh.getCoord(),this.id);
				if(occupied == false){
					obj.add(neigh);
				}
			}
		}

		if(obj.size()<2)obj=getVisibleFoodPatches(2*Parameter.visual_range);

		return obj;
	}

	private ArrayList<Cell> getVisibleFoodPatches(double f, ArrayList<Cell> food){

		ArrayList<Cell> obj = new ArrayList<Cell>();
		for(Cell c: food){
			if(c.getCoord().distance(this.getCoord())<f){
				boolean occupied = SimUtils.occupied(c.getCoord(),this.id);
				if(occupied == false){
					obj.add(c);
				}
			}
		}

		if(obj.size()<2)obj=getVisibleFoodPatches(2*Parameter.visual_range, food);

		return obj;
	}

	private Cell getCurrentFoodPatch(){

		Iterable<Cell> objectsInArea = null;
		Envelope envelope = new Envelope();
		envelope.init(this.getCoord());
		envelope.expandBy(Math.max(Parameter.foodBuffer,4));
		objectsInArea = ModelSetup.getGeog().getObjectsWithin(envelope,Cell.class);
		envelope.setToNull();

		Cell myCell = null;
		double minDist = 999999;
		while(objectsInArea.iterator().hasNext()){
			Cell neigh = objectsInArea.iterator().next();
			if(neigh.getCoord().distance(this.getCoord())<=Parameter.foodBuffer){
				double dist = neigh.getCoord().distance(this.getCoord());
				if(dist<minDist){
					myCell = neigh;
					minDist=dist;
				}	
			}

		}

		return myCell;
	}

	private synchronized void move(RealVector direction,boolean isCellSite){
		//this.setFacing(direction);
		MoveUtils.moveTo((Primate)this, direction);
		this.setCoord(ModelSetup.getAgentGeometry(this).getCoordinates()[0]);
	}


}