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
	static List<Baboon> primateList;
	ArrayList<Double> speeds,spreads,ratio,perpen,shapeRatio;
	ArrayList<Coordinate> centers;
	static ArrayList<ArrayList> indPositions,foodQuant;
	static ArrayList<Double> foodPosX,foodPosY;
	static ArrayList<Double> centerDist_path,indDist_path,foragingIntake_mean,foragingIntake_sd;
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
		centerDist_path = new ArrayList<Double>();
		indDist_path = new ArrayList<Double>();
		foragingIntake_mean = new ArrayList<Double>();
		foragingIntake_sd = new ArrayList<Double>();


		Collections.sort(primateList,new CustomComparator());
		for(Primate p : primateList){
			//System.out.println("id "+ p.getId());
		}

		for(Primate p : primateList){
			//System.out.println("id "+ p.getId());
		}

		//creating a file to store the output of the counts
		try {
			summaryStats_out = new BufferedWriter(new FileWriter("data/summary_stats.csv",true));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	/********************************	
	 * 								*
	 *         Stepper				*
	 * 								*
	 *******************************/

	@ScheduledMethod(start=26, interval = 1,priority=0)
	public void step(){

		//	to do every tick

		

		//	to do every minute
		if(RunEnvironment.getInstance().getCurrentSchedule().getTickCount()%(Parameter.recordingFreq)==0){

			getDistanceFromPath();
			getForagingEfficiency();
			
		}

		//	to do every 1000 steps
		if(RunEnvironment.getInstance().getCurrentSchedule().getTickCount()%(1000)==0){
			//hour++;
		}

		//	to do at the end of the simulation
		if(RunEnvironment.getInstance().getCurrentSchedule().getTickCount()>=Parameter.stepsPerDay*Parameter.endDay){
			executor.shutdown();
			RunEnvironment.getInstance().endAt(this.runEnvironment.getCurrentSchedule().getTickCount());
			recordSummaryStats();
			
			System.out.println("End of Sim");
			
		}
	}


	/****************************************Methods****************************************************************/

	private static void getDistanceFromPath(){

		//distance based on group center
		Coordinate groupCenter = calculateGroupCenter();
		centerDist_path.add(getDistanceBetweenCoordAndPath(groupCenter));


		//distance based on average individual distance
		double meanDist=0;
		int count=0;
		for(Coordinate cc : getIndividualCoordinates()){
			meanDist = meanDist + getDistanceBetweenCoordAndPath(cc);
			count++;
		}
		meanDist = meanDist/(double)count;
		indDist_path.add(meanDist);
	}

	private static double getDistanceBetweenCoordAndPath(Coordinate cc){

		double minDistance = 9999;
		for(markerPoint pathCell : ModelSetup.pathCells){
			double dist = pathCell.getCoord().distance(cc);
			if(dist<minDistance)minDistance=dist;
		}

		return minDistance;
	}

	private static void getForagingEfficiency(){

		//get individual mean and sd of foraging efficiency
		double meanF=0,sdF=0;
		for(Baboon bb : ModelSetup.orderedP){
			meanF = meanF +bb.getEnergyIntake();
		}
		meanF = meanF / (double)ModelSetup.orderedP.size();

		for(Baboon bb : ModelSetup.orderedP){
			sdF = sdF + Math.pow(meanF-bb.getEnergyIntake(), 2);
		}
		sdF = sdF / ((double)ModelSetup.orderedP.size()-1);

		foragingIntake_mean.add(meanF);
		foragingIntake_sd.add(sdF);
	}



	private static Coordinate calculateGroupCenter(){
		//find the center of this group
		double xCoordAvg=0;
		double yCoordAvg=0;
		double numbPrimates=0;

		//get all individuals from focal group
		ArrayList<Primate> group = new ArrayList<Primate>();
		for(Primate p : primateList){
			group.add(p);
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

	private static ArrayList<Coordinate> getIndividualCoordinates(){
		ArrayList<Coordinate> positions = new ArrayList<Coordinate>();
		for (Primate p: ModelSetup.orderedP){
			positions.add(p.getCoord());
		}
		return positions;
	}

	private double calculateSpread(Coordinate groupCenter){

		double MSE_dist=-1;
		int numbPrimates=0;

		ArrayList<Primate> group = new ArrayList<Primate>();
		for(Primate p : primateList){
			group.add(p);
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

		double avg_centerDistance=0,avg_individualDistance=0,f_intake=0,fsd_intake=0;

		//distance from center to path (mean)
		for(Double d:centerDist_path){
			avg_centerDistance = avg_centerDistance + d;
		}
		avg_centerDistance = avg_centerDistance/(double)centerDist_path.size();

		//distance from individuals to path (mean)
		for(Double d:indDist_path){
			avg_individualDistance = avg_individualDistance + d;
		}
		avg_individualDistance = avg_individualDistance/(double)indDist_path.size();

		//Mean foraging intake from all individuals (mean)
		for(Double d:foragingIntake_mean){
			f_intake = f_intake + d;
		}
		f_intake = f_intake/(double)foragingIntake_mean.size();

		//Variation in foraging intake from all individuals (mean)
		for(Double d:foragingIntake_sd){
			fsd_intake = fsd_intake + d;
		}
		fsd_intake = fsd_intake/(double)foragingIntake_sd.size();


		try {
			//setup headers
			//summaryStats_out.append("centerPathDist,indPathDist,ForagingIntake,foraginIntakeSD");
			//summaryStats_out.newLine();

			//write out data
			summaryStats_out.append(((Double)avg_centerDistance).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Double)avg_individualDistance).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Double)f_intake).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Double)fsd_intake).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Integer)Parameter.groupSize).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Integer)Parameter.influenceType).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Integer)Parameter.addPath).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Double)Parameter.corePer).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Double)Parameter.pathFood).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Integer)Parameter.pathWidth).toString());
			summaryStats_out.append(",");
			summaryStats_out.append(((Double)Parameter.stepsPerDay).toString());
			summaryStats_out.append(",");
			
		//summaryStats_out.newLine();
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
