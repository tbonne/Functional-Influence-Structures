package tools;

import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.linear.RealVector;
import org.geotools.referencing.GeodeticCalculator;

import repast.simphony.random.RandomHelper;
import LHP.Baboon;
import LHP.Cell;
import LHP.Parameter;
import LHP.Primate;
import LHP.ModelSetup;

import com.vividsolutions.jts.geom.Coordinate;

public final strictfp  class MoveUtils {

	/*private synchronized static boolean moveToS(Primate primate, Coordinate c,int attemptCount){

		if(primate.getId()==8){
			System.out.println("");
		}

		boolean retval=false;
		Coordinate newCoord = c;
		int count = attemptCount;

		double angle=getAngle(primate.getCoord(),newCoord);
		
		//check to make sure not moving more than the max distance
		double moveDist = (primate.getCoord().distance(newCoord));
		if (moveDist>Parameter.maxDistancePerStep)moveDist=Parameter.maxDistancePerStep;

		//check to make sure nan is represented as 0 distance
		if (Double.isNaN(moveDist) || moveDist==0.0){
			//do not move, 0 distance
		} else {
			//if my way is not blocked 
			if(primate.blocked()==false){

				//am i facing my target?
				if(Math.abs(angle-primate.getFacing())<=1){

					//target coordinate for this step
					Coordinate targetCoord = new Coordinate(primate.getCoord().x+Math.cos(Math.toRadians(angle))*moveDist,primate.getCoord().y+Math.sin(Math.toRadians(angle))*moveDist);

					//check to make sure not to move onto another individual
					boolean occupied=checkIfWayIsBlocked(primate,targetCoord);
					//boolean overlap=checkForOverlap(primate);

					//The way is clear
					if(occupied == false){
						//Not blocked, perform move
						primate.setFacing(angle);
						ModelSetup.getGeog().moveByDisplacement(primate, Math.cos(Math.toRadians(angle))*moveDist, Math.sin(Math.toRadians(angle))*moveDist);
						//primates make noise when moving
						if(moveDist>0)primate.setNoise(1);

						//the way is blocked
					} else if (occupied==true){
						//look for a way around
						int i=0;
						if(RandomHelper.nextDouble()>0.5){
							i=-1;
						} else {
							i=1;
						}
						primate.setFacing(primate.getFacing()+(10*i));
						primate.setBlocked(true);
					}
					//not facing the right direction, turn to face target
				} else {
					primate.setFacing(AngleUtils.turnToFaceTarget(primate.getFacing(),angle));
					ModelSetup.getGeog().moveByDisplacement(primate, Math.cos(Math.toRadians(primate.getFacing()))*moveDist, Math.sin(Math.toRadians(primate.getFacing()))*moveDist);
					if(moveDist>0)primate.setNoise(1);
					//primate.setFacing(angle);
				}
			} else if (primate.blocked() == true){
				primate.setFacing(angle);
				Coordinate targetCoord = new Coordinate(primate.getCoord().x+Math.cos(Math.toRadians(primate.getFacing()))*Parameter.maxDistancePerStep,primate.getCoord().y+Math.sin(Math.toRadians(primate.getFacing()))*Parameter.maxDistancePerStep);
				boolean occupied=checkIfWayIsBlocked(primate,targetCoord);

				if(occupied==true){
					double turn = RandomHelper.nextDoubleFromTo(-10,10);
					primate.setFacing(primate.getFacing()+turn);
					ModelSetup.getGeog().moveByDisplacement(primate, Math.cos(Math.toRadians(primate.getFacing()))*Parameter.maxDistancePerStep, Math.sin(Math.toRadians(primate.getFacing()))*Parameter.maxDistancePerStep);
					if (moveDist>0)retval=true;
				} else {
					ModelSetup.getGeog().moveByDisplacement(primate, Math.cos(Math.toRadians(primate.getFacing()))*Parameter.maxDistancePerStep, Math.sin(Math.toRadians(primate.getFacing()))*Parameter.maxDistancePerStep);
					primate.setBlocked(false);
				}
			} 
		}
		
		return retval;
	}*/


	/*public synchronized static boolean moveTo(Primate primate, Coordinate c){

		boolean retval=false;
		Coordinate newCoord = c;

		double angle=getAngle(primate.getCoord(),newCoord);
		
		double moveDist = (primate.getCoord().distance(newCoord));
		
		//check to make sure not moving more than the max distance
		if (moveDist>Parameter.maxDistancePerStep)moveDist=Parameter.maxDistancePerStep;

		//Move in the direction of target
		if(moveDist>0 &&Double.isNaN(moveDist)==false){
			primate.setFacing(angle);
			ModelSetup.getGeog().moveByDisplacement(primate, Math.cos(Math.toRadians(angle))*moveDist, Math.sin(Math.toRadians(angle))*moveDist);
			retval=true;
		}
				
		return retval;
	}*/
	
	public synchronized static boolean moveTo(Primate primate, RealVector vecDest){

		boolean retval=false;

		ModelSetup.getGeog().moveByDisplacement(primate, vecDest.getEntry(0), vecDest.getEntry(1));
		//ModelSetup.getGeog().moveByDisplacement(primate, 0, -0.1);
				
		return retval;
	}
	
	/*
	public static boolean checkIfWayIsBlocked(Primate primate,Coordinate targetCoord){
		boolean occupied=false;
		for(Primate p:ModelSetup.primatesAll){
			if(p.getId()!=primate.getId()){
				double neighDist = p.getCoord().distance(targetCoord);
				if(neighDist<=Parameter.bodyRadius){
					occupied = true;
					//System.out.println("occupied!");
				}
			}
		}
		return occupied;
	}
	
	public static Coordinate checkIfWayIsBlockedC(Primate primate,Coordinate targetCoord){
		Coordinate block=null;
		for(Primate p:ModelSetup.primatesAll){
			if(p.getId()!=primate.getId()){
				double neighDist = p.getCoord().distance(targetCoord);
				if(neighDist<=Parameter.bodyRadius){
					block=p.getCoord();
					//System.out.println("occupied!");
				}
			}
		}
		return block;
	}

	public static boolean checkForOverlap(Primate primate){
		boolean overlap = false;
		for(Primate p:ModelSetup.primatesAll){
			if(p.getId()!=primate.getId()){
				double neighDist = p.getCoord().distance(primate.getCoord());
				if(neighDist<=Parameter.bodyRadius){
					overlap = true;
					//System.out.println("overlap!");
				}
			}
		}
		return overlap;
	}*/


	public static double getAngle(Coordinate source, Coordinate target) {
		float angle = (float) (Math.atan2(target.y - source.y, target.x - source.x));

		if(angle < 0){
			angle += 2*Math.PI;
		}
		return Math.toDegrees(angle);
	}
	
	public static double getAngle(Coordinate source, Coordinate target, boolean rads) {
		float angle = (float) (Math.atan2(target.y - source.y, target.x - source.x));

		if(angle < 0){
			angle += 2*Math.PI;
		}
		return angle;
	}


}
