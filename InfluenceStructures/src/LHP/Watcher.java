package LHP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.ShapefileWriter;

public class Watcher {

	Executor executor;
	RunEnvironment runEnvironment;
	BufferedWriter summaryStats_out,individualMovements_out,foodQuant_out,foodPos_out;//locations
	List<Baboon> primateList;
	ArrayList<Double> speeds,spreads,ratio,perpen,shapeRatio;
	ArrayList<Coordinate> centers;
	static ArrayList<ArrayList> indPositions,foodQuant;
	static ArrayList<Double> foodPosX,foodPosY;
	double hour;


	//output stats
	Coordinate lastCenter;

	public Watcher(Executor exe){

		executor = exe;
		runEnvironment = RunEnvironment.getInstance();
		primateList = ModelSetup.primatesAll;
		speeds = new ArrayList<Double>();
		ratio = new ArrayList<Double>();
		spreads = new ArrayList<Double>();
		centers = new ArrayList<Coordinate>();
		indPositions = new ArrayList<ArrayList>();
		hour = 1;
		foodPosX = new ArrayList<Double>();
		foodPosY = new ArrayList<Double>();
		foodQuant = new ArrayList<ArrayList>();
		perpen = new ArrayList<Double>();
		shapeRatio = new ArrayList<Double>();


		Collections.sort(primateList,new CustomComparator());
		for(Primate p : primateList){
			//System.out.println("id "+ p.getId());
		}

		for(Primate p : primateList){
			//System.out.println("id "+ p.getId());
		}

		//creating a file to store the output of the counts
		try {
			//locations = new BufferedWriter(new FileWriter("C:/Users/t-work/Dropbox/Projects_with_LouisePeter/Project_Babbon_group_cordination_models_(specificFocus)/Project_simple_to_complex_models_min_requirments/runingABM/Primate_locations"+Parameter.decisionMaking+"_"+Parameter.associationSize+".csv",false));
			summaryStats_out = new BufferedWriter(new FileWriter("C:/Users/t-work/Dropbox/Tyler_MSS/GeladaOMU_Model/results/summary_stats_"+Parameter.numbOfGroups+"_OMU.csv",false));
			//individualMovements_out = new BufferedWriter(new FileWriter("C:/Users/t-work/Dropbox/Projects_with_LouisePeter/Project_Babbon_group_cordination_models_(specificFocus)/Project_simple_to_complex_models_min_requirments/runingABM/movement_stats.csv",false));
			//foodQuant_out = new BufferedWriter(new FileWriter("C:/Users/t-work/Dropbox/Projects_with_LouisePeter/Project_Babbon_group_cordination_models_(specificFocus)/Project_simple_to_complex_models_min_requirments/runingABM/food_quant.csv",false));
			//foodPos_out = new BufferedWriter(new FileWriter("C:/Users/t-work/Dropbox/Projects_with_LouisePeter/Project_Babbon_group_cordination_models_(specificFocus)/Project_simple_to_complex_models_min_requirments/runingABM/food_pos.csv",false));
		} catch (IOException e) {
			e.printStackTrace();
		}

		//logNDist = new LogNormalDistribution(scale,shape);
	}


	/********************************	
	 * 								*
	 *         Stepper				*
	 * 								*
	 *******************************/

	@ScheduledMethod(start=26, interval = 1,priority=0)
	public void step(){

		//	to do every tick

		//countSteps=countSteps+1;
		//if(revisit==null)revisit=Math.round(logNDist.sample()); can't down-sample this way as it assumes that all positions are recorded at the same time... no interpolation error introduced.
		//if(countSteps==revisit){
		//recordIndividualPositions();
		//recordFood();
		//	revisit=null;
		//	countSteps=0;
		//}


		//	to do every minute
		if(RunEnvironment.getInstance().getCurrentSchedule().getTickCount()%(Parameter.recordingFreq)==0){
			Coordinate groupCenter = calculateGroupCenter();
			centers.add(groupCenter);
			speeds.add(calculateSpeed(groupCenter));
			spreads.add(calculateSpread(groupCenter));
			//ratio.add(calculateLegthWidthRatio());
			perpen.add(this.calculateDiffAngle(this.getDirectionOfTravel().getEntry(1)/this.getDirectionOfTravel().getEntry(0), this.calculateSlopeOfEllipse().getEntry(1)/this.calculateSlopeOfEllipse().getEntry(0)));
			shapeRatio.add(this.calculateShape(Math.atan2(this.calculateSlopeOfEllipse().getEntry(1), this.calculateSlopeOfEllipse().getEntry(0))));
			System.out.println("step");
		}

		//	to do every 1000 steps
		if(RunEnvironment.getInstance().getCurrentSchedule().getTickCount()%(1000)==0){
			//System.out.println("End of hour "+ hour);
			//System.out.println(" ");
			hour++;
		}

		//	to do at the end of the simulation
		if(RunEnvironment.getInstance().getCurrentSchedule().getTickCount()>=Parameter.stepsPerDay*Parameter.endDay){
			executor.shutdown();
			RunEnvironment.getInstance().endAt(this.runEnvironment.getCurrentSchedule().getTickCount());
				recordSummaryStats();
			//	recordIndPosisions();
			//	recordFoodPo();
			//	recordFood_txt();
			//	closeBufferWiter();
			//	closeBufferWiter2();
			//	printSocialCounts();
			System.out.println("End of day");
			//exportLandscape();
		}
	}

	//private void printSocialCounts(){
	//	
	//	for(Primate p : ModelSetup.orderedP){
	//		System.out.println(p.getId()+" "+p.socialCount);
	//	}
	//	
	//}

	private double calculateSpeed(Coordinate groupCenter){

		double speed=-1;

		if(lastCenter!=null){
			speed = groupCenter.distance(lastCenter)/1;
		} 

		lastCenter = groupCenter;

		return speed;
	}

	

	private RealVector getDirectionOfTravel(){

		double direction=0;
		//get all individuals from focal group
		ArrayList<Primate> group = new ArrayList<Primate>();
		for(Primate p : primateList){
			if(p.getMyGroup()==0)group.add(p);
		}

		//get direction of travel
		double xsum=0, ysum=0;
		for(Primate g: group){
			xsum=xsum+g.getFacing().getEntry(0);
			ysum=ysum+g.getFacing().getEntry(1);
		}
		RealVector directionOfTravel = new ArrayRealVector(2,0);
		directionOfTravel.addToEntry(0, xsum/(double)group.size());
		directionOfTravel.addToEntry(1, ysum/(double)group.size());

		//System.out.println("Direction of travel of the group = "+directionOfTravel.toString());
		directionOfTravel.unitize();
		System.out.println("Direction of travel of the group unit vector = "+directionOfTravel.toString());
		//direction = Math.atan2(directionOfTravel.getEntry(1), directionOfTravel.getEntry(0));
		
		
		return directionOfTravel;
	}

	private RealVector calculateSlopeOfEllipse(){

		//get all individuals from focal group
		ArrayList<Primate> group = new ArrayList<Primate>();
		for(Primate p : primateList){
			if(p.getMyGroup()==0)group.add(p);
		}

		//getGroupCenter
		Coordinate center = calculateGroupCenter();

		//calculate angle of ellipse (source arcgis standard deviation ellipse)
		double A = 0,B=0,C=0;

		double xdist2 = 0,ydist2=0, ydistxdist=0;
		for(Primate g: group){
			//individual location based on the center of the group
			xdist2 = xdist2+ Math.pow(g.getCoord().x-center.x,2);
			ydist2 = ydist2+ Math.pow(g.getCoord().y-center.y,2);
			ydistxdist = ydistxdist + (g.getCoord().x-center.x)*(g.getCoord().y-center.y);
		}

		A = xdist2-ydist2;

		B = Math.pow( Math.pow(A,2) + 4*(Math.pow(ydistxdist, 2)), 0.5);

		C = 2*ydistxdist;

		double slope = Math.atan((A+B)/C);
		
		RealVector slopeE = new ArrayRealVector(2,0);
		slopeE.addToEntry(1, A+B);
		slopeE.addToEntry(0, C); 

		return slopeE;
	}
	
	private double calculateDiffAngle(double slopeD, double slopeE){
		double angleD=Math.atan(slopeD);
		//slopeD = 10000;
		double angleE=Math.atan(slopeE);
		double diff = angleD-angleE;
		double diff2 = Math.atan((slopeD-slopeE)/(1+slopeD*slopeE));
		diff2 =  ((Math.PI/2.0)-Math.abs(diff2))/(Math.PI/2.0);
		if(Double.isNaN(diff2)==true && diff == 0)diff2=0;
		
		System.out.println("SlopeD = "+slopeD+",  SlopeE = "+slopeE);
		System.out.println("perpend = "+diff2);
		
		return diff2;
	}

	private double calculateShape(Double angleE){

		//get all individuals from focal group
		ArrayList<Primate> group = new ArrayList<Primate>();
		for(Primate p : primateList){
			if(p.getMyGroup()==0)group.add(p);
		}

		//getGroupCenter
		Coordinate center = calculateGroupCenter();

		//calculate sd of ellipse (source arcgis standard deviation ellipse)
		double xd = 0,yd=0;
		for(Primate g: group){
			xd = xd + Math.pow( (g.getCoord().x-center.x)*Math.cos(angleE)-(g.getCoord().y-center.y)*Math.sin(angleE),2);
			yd = yd + Math.pow( (g.getCoord().x-center.x)*Math.sin(angleE)-(g.getCoord().y-center.y)*Math.cos(angleE),2);
		}
		
		xd = Math.pow(xd/(double)group.size(),0.5);
		yd = Math.pow(yd/(double)group.size(),0.5);
		
		double sl = 0;
		if(xd>=yd){
			sl=yd/xd;
		} else{
			sl=xd/yd;
		}
		
		System.out.println("shape = "+(1-sl));
		
		return 1-sl;

	}

	/*//Calculate the ratio of length to width based on max
	double maxX=-99999, minX=99999, maxY = -999999, minY=999999;
	for(RealVector v : rotatedPoints){

		double x = v.getEntry(0);
		double y = v.getEntry(1);

		if(x>maxX)maxX=x;
		if(x<minX)minX=x;
		if(y>maxY)maxY=y;
		if(y<minY)minY=y;

	}

	lengthWidth = (maxY-minY)/(maxX-minX);
	
	private double calculateLegthWidthRatio(){

		double lengthWidth = -1;

		//get all individuals from focal group
		ArrayList<Primate> group = new ArrayList<Primate>();
		for(Primate p : primateList){
			if(p.getMyGroup()==0)group.add(p);
		}

		//getGroupCenter
		Coordinate center = calculateGroupCenter();

		//get direction of travel
		double xsum=0, ysum=0;
		for(Primate g: group){
			xsum=xsum+g.getFacing().getEntry(0);
			ysum=ysum+g.getFacing().getEntry(1);
		}
		RealVector directionOfTravel = new ArrayRealVector(2,0);
		directionOfTravel.addToEntry(0, xsum/(double)group.size());
		directionOfTravel.addToEntry(1, ysum/(double)group.size());

		System.out.println("Direction of travel of the group = "+directionOfTravel.toString());
		directionOfTravel.unitize();
		System.out.println("Direction of travel of the group unit vector = "+directionOfTravel.toString());

		//subtract all the points towards the center and rotate
		ArrayList<RealVector> rotatedPoints = new ArrayList<RealVector>();
		for(Primate g: group){

			//individual location based on the center of the group
			double x = g.getCoord().x-center.x;
			double y = g.getCoord().y-center.y;

			//calculate angle of rotation
			double fie = Math.atan2(y, x);
			fie = -Math.atan2(directionOfTravel.getEntry(1),directionOfTravel.getEntry(0));

			//Rotate point
			RealVector p1 = new ArrayRealVector(2,0);
			p1.addToEntry(0, x*Math.cos(fie) - y*Math.sin(fie));
			p1.addToEntry(1, y*Math.cos(fie) + x*Math.sin(fie));
			rotatedPoints.add(p1);

		}

		//Calculate the ratio of length to width based on SD
		double SDx=0,SDy=0;//mean = 0;
		for(RealVector v : rotatedPoints){
			SDx = SDx+Math.pow(0-v.getEntry(0),2);
			SDy = SDy+Math.pow(0-v.getEntry(1),2);
		}

		if(SDy>=SDx){
			lengthWidth = ( Math.pow(SDy/(double)rotatedPoints.size(),0.5) )/( Math.pow(SDx/(double)rotatedPoints.size(),0.5) );//N cancels out
		}else{
			lengthWidth = ( Math.pow(SDx/(double)rotatedPoints.size(),0.5) )/( Math.pow(SDy/(double)rotatedPoints.size(),0.5) );//N cancels out
		}
		
		return lengthWidth;

	}
	 */

	private static void recordFood(){
		ArrayList foodT = new ArrayList<Double>();
		for(Cell c: ModelSetup.allCells){
			foodT.add(c.resources);
		}
		foodQuant.add(foodT);
	}

	private static void recordFoodPo(){

		for(Cell c: ModelSetup.allCells){
			foodPosX.add(c.getCoord().x);
			foodPosY.add(c.getCoord().y);
		}
	}

	private void recordFood_txt(){
		try {
			int count=0;
			for(int i = 0 ; i<foodQuant.size();i++){
				ArrayList fq = foodQuant.get(i);
				for(int j = 0 ; j < fq.size() ; j++ ){
					foodQuant_out.append(((Double)fq.get(j)).toString());
					foodQuant_out.append(",");
				}
				foodQuant_out.newLine();
			}

			for(int k = 0 ; k < foodPosX.size();k++){
				foodPos_out.append(((Double)foodPosX.get(k)).toString());
				foodPos_out.append(",");
				foodPos_out.append(((Double)foodPosY.get(k)).toString());
				foodPos_out.newLine();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		try {

			foodQuant_out.flush();
			foodQuant_out.close();

			foodPos_out.flush();
			foodPos_out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private Coordinate calculateGroupCenter(){
		//find the center of this group
		double xCoordAvg=0;
		double yCoordAvg=0;
		double numbPrimates=0;

		//get all individuals from focal group
		ArrayList<Primate> group = new ArrayList<Primate>();
		for(Primate p : primateList){
			if(p.getMyGroup()==0)group.add(p);
		}

		//calculate their average x and y coordinates
		for (Primate p: group){
			xCoordAvg = xCoordAvg + p.getCoord().x;
			yCoordAvg = yCoordAvg + p.getCoord().y;
			numbPrimates++;
		}

		xCoordAvg = xCoordAvg/numbPrimates;
		yCoordAvg = yCoordAvg/numbPrimates;

		//return the center coordinate
		return new Coordinate(xCoordAvg,yCoordAvg);
	}

	private static void recordIndividualPositions(){
		ArrayList<Double> positions = new ArrayList<Double>();
		for (Primate p: ModelSetup.orderedP){
			positions.add(p.getCoord().x);
			positions.add(p.getCoord().y);
		}
		indPositions.add(positions);
	}

	private double calculateSpread(Coordinate groupCenter){

		double MSE_dist=-1;
		int numbPrimates=0;
		
		ArrayList<Primate> group = new ArrayList<Primate>();
		for(Primate p : primateList){
			if(p.getMyGroup()==0)group.add(p);
		}

		for (Primate p: group){
			MSE_dist = MSE_dist + Math.pow(p.getCoord().distance(groupCenter),2);
			numbPrimates++;
		}

		return MSE_dist/numbPrimates;
	}

	private void recordIndPosisions(){
		try {
			int count=0;
			for(int i = 0 ; i<indPositions.size();i++){
				ArrayList<Double> pos = indPositions.get(i);
				for(int j = 0 ; j<pos.size();j++){
					//x
					individualMovements_out.append(((Double)pos.get(j)).toString());
					individualMovements_out.append(",");
					j++;
					//y
					individualMovements_out.append(((Double)pos.get(j)).toString());
					individualMovements_out.append(",");
					//id
					individualMovements_out.append(((Integer)count).toString());
					individualMovements_out.append(",");
					count++;
					//time
					individualMovements_out.append(((Integer)i).toString());

					//new line
					individualMovements_out.newLine();
				}
				count=0;
				//individualMovements_out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		flushBufferWriter2();
	}



	private void recordSummaryStats() {

		double avg_speed=0,avg_spread=0;

		for(Double d:spreads){
			avg_spread = avg_spread + d;
		}
		avg_spread = avg_spread/spreads.size();
		for(Double d:speeds){
			avg_speed = avg_speed + d;
		}
		avg_speed = avg_speed / speeds.size();

		System.out.println("speed = "+avg_speed+"  spread = "+avg_spread);
		
		

		try {
			//setup headers
			summaryStats_out.append("time,x,y,speed,spread,perpendicularity,shape,rankF");
			summaryStats_out.newLine();
			
			//write out data
			for(int i = 0 ; i<centers.size();i++){
				summaryStats_out.append(((Double)(i*Parameter.recordingFreq)).toString());
				summaryStats_out.append(",");
				summaryStats_out.append(((Double)centers.get(i).x).toString());
				summaryStats_out.append(",");
				summaryStats_out.append(((Double)centers.get(i).y).toString());
				summaryStats_out.append(",");
				summaryStats_out.append(((Double)speeds.get(i)).toString());
				summaryStats_out.append(",");
				summaryStats_out.append(((Double)spreads.get(i)).toString());
				summaryStats_out.append(",");
				summaryStats_out.append(((Double)this.perpen.get(i)).toString());
				summaryStats_out.append(",");
				summaryStats_out.append(((Double)this.shapeRatio.get(i)).toString());
				summaryStats_out.append(",");
				summaryStats_out.append(  ((Double)((this.shapeRatio.get(i))*(this.perpen.get(i)))).toString());
				summaryStats_out.append(",");
				
				summaryStats_out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		flushBufferWriter();
	}


	private void flushBufferWriter(){
		try {
			summaryStats_out.newLine();
			summaryStats_out.flush();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void flushBufferWriter2(){
		try {
			individualMovements_out.newLine();
			individualMovements_out.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void closeBufferWiter(){
		try {
			summaryStats_out.flush();
			summaryStats_out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void closeBufferWiter2(){
		try {
			individualMovements_out.flush();
			individualMovements_out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public class CustomComparator implements Comparator<Primate> {
		@Override
		public int compare(Primate o1, Primate o2) {
			int retval = 1;
			if(o1.getId()<o2.getId())retval=0;
			return retval;
		}
	}

	private void exportLandscape(){
		System.out.println("Exporting landscape to be analyzed");

		File point = new File("C:/Users/t-work/Desktop/GIS/BehaviourExtraction/Landscape"+".shp");

		ShapefileWriter shapeOut = new ShapefileWriter(ModelSetup.getGeog());
		try {
			shapeOut.write(ModelSetup.getGeog().getLayer(Cell.class).getName(), point.toURI().toURL());
		} catch (MalformedURLException e ) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NullPointerException n){
			n.printStackTrace();
		}

	}

}
