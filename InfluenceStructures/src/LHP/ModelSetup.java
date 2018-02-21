package LHP;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import org.geotools.referencing.GeodeticCalculator;

import jsc.distributions.Lognormal;
import jsc.distributions.Normal;
import cern.jet.random.Beta;
import cern.jet.random.Uniform;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.context.space.gis.GeographyFactory;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import tools.NetworkUtils;
import tools.SimUtils;

//This file builds the model: creating the environment and populating it with agents
public class ModelSetup implements ContextBuilder<Object>{

	private static Context mainContext;
	private static Geography geog;
	private static int sitesAdded;
	private static int previousListSize;
	private static double resAdded;
	public static int primatesAdded;
	public static ArrayList<Baboon> primatesAll,orderedP;
	public static Iterable<Primate> primateIterator;
	public static ArrayList<Cell> cellsToUpdate;
	public static ArrayList<Cell> removeCellsToUpdate;
	public static ArrayList<Cell> allCells;
	public static ArrayList<markerPoint> pathCells;
	public static double timeRecord;
	public static double timeRecord_start;
	public static GeodeticCalculator gc;

	public Context<Object> build(Context<Object> context){

		System.out.println("Running movement pattern simulation model: group level foraging efficency");

		/********************************
		 * 								*
		 * initialize model parameters	*
		 * 								*
		 * ******************************/

		sitesAdded = 0;
		resAdded=0;
		primatesAdded = 0;
		primatesAll= new ArrayList<Baboon>();
		orderedP = new ArrayList<Baboon>();
		mainContext = null;
		Parameter parameter = new Parameter();
		geog=null;
		cellsToUpdate = new ArrayList<Cell>();
		removeCellsToUpdate = new ArrayList<Cell>();
		allCells = new ArrayList<Cell>();
		pathCells = new ArrayList<markerPoint>();
		mainContext = context; //static link to context
		timeRecord = System.currentTimeMillis();
		timeRecord_start = System.currentTimeMillis();

		//get parameters from csv file
		if(parameter.getParamsCSV)resetP();

		/****************************
		 * 							*
		 * Building the landscape	*
		 * 							*
		 * *************************/

		//Create Geometry factory; used to create gis shapes (points=primates; polygons=resources)
		GeometryFactory fac = new GeometryFactory();


		//x and y dims of the map file
		int xdim = Parameter.landscapeWidth;
		int ydim = Parameter.landscapeHeight;

		//Create Geography/GIS 
		GeographyParameters<Object> params= new GeographyParameters<Object>();
		GeographyFactory factory = GeographyFactoryFinder.createGeographyFactory(null);
		geog = factory.createGeography(Parameter.geog, context, params);
		geog.setCRS("EPSG:32636"); //WGS 84 / UTM zone 36N EPSG:32636
		gc = new GeodeticCalculator(geog.getCRS());

		//Add Resources to the environment ****************************************

		System.out.println("adding resources to the environment");
		Beta beta = new Beta(Parameter.envHomogen,Parameter.envHomogen, Parameter.mt);
		Uniform xDist = RandomHelper.createUniform(0, xdim);
		Uniform yDist = RandomHelper.createUniform(0, ydim);
		int count=0;
		for(int i=0;i<Parameter.foodDensity*xdim*ydim;i++){
			Cell cell = new Cell(context,xDist.nextDouble(),yDist.nextDouble(),beta.nextDouble(),count++);
			sitesAdded++;
			allCells.add(cell);
		}

		//increase area of the foraging trail (create 200m buffer around original surface) 
		System.out.println("adding more resources to the environment");
		Uniform xDist_buff_top = RandomHelper.createUniform(0, xdim+200); //top buffer area
		Uniform yDist_buff_top = RandomHelper.createUniform(ydim, ydim+200);
		for(int i=0;i<Parameter.foodDensity*(xdim+200)*200;i++){
			Cell cell = new Cell(context,xDist_buff_top.nextDouble(),yDist_buff_top.nextDouble(),beta.nextDouble(),count++);
			sitesAdded++;
			allCells.add(cell);
		}
		Uniform xDist_buff_right = RandomHelper.createUniform(xdim, xdim+200); //right buffer area
		Uniform yDist_buff_right = RandomHelper.createUniform(0,ydim);
		for(int i=0;i<Parameter.foodDensity*200*ydim;i++){
			Cell cell = new Cell(context,xDist_buff_right.nextDouble(),yDist_buff_right.nextDouble(),beta.nextDouble(),count++);
			sitesAdded++;
			allCells.add(cell);
		}
		Uniform xDist_buff_bottom = RandomHelper.createUniform(0, xdim+200); //bottom buffer area
		Uniform yDist_buff_bottom = RandomHelper.createUniform(0,-200);
		for(int i=0;i<Parameter.foodDensity*(xdim+200)*200;i++){
			Cell cell = new Cell(context,xDist_buff_bottom.nextDouble(),yDist_buff_bottom.nextDouble(),beta.nextDouble(),count++);
			sitesAdded++;
			allCells.add(cell);
		}


		//Add high density paths to the environment ****************************************

		if(Parameter.addPath == 1){

			//add Resources to the environment (high density path)
			Uniform yDist2 = RandomHelper.createUniform(0, ydim);

			cern.jet.random.Normal err = RandomHelper.createNormal(0, Parameter.pathWidth);
			for(int i=0;i<Parameter.pathFood*xdim*ydim;i++){
				double ySample = yDist2.nextDouble();
				double xSample = xdim+1;
				while (xSample > xdim || xSample < 0){
					xSample = Math.pow( ySample, 0.25) ;
					xSample = xSample / Math.pow(ydim,0.25);
					xSample = xSample * xdim ;
				}

				Cell cell = new Cell(context,xSample+ err.nextDouble(),ySample+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
			}

			//add marker points with no error (act as line string for distance calculations)
			System.out.println("adding marker points to the environment"); 
			for(int i=0;i<ydim;i++){
				double ySample = i;//yDist2.nextDouble();
				double xSample = xdim+1;
				while (xSample > xdim || xSample < 0){
					//ySample = -1*Math.pow((xSample-xdim/2)/2,2)+ydim/2+err.nextDouble();
					xSample = Math.pow( ySample, 0.25);
					xSample = xSample / Math.pow(ydim,0.25);
					xSample = xSample * xdim ;
				}

				markerPoint mp = new markerPoint(context,xSample,ySample);
				pathCells.add(mp);
			}

		} else if (Parameter.addPath==2){
			
			//path length
			double pathLength = ydim/(2*2);
			
			cern.jet.random.Normal err = RandomHelper.createNormal(0, Parameter.pathWidth);
			cern.jet.random.Uniform whichPath = RandomHelper.createUniform(1,2);
			cern.jet.random.Uniform wherePath = RandomHelper.createUniform(0,pathLength);
			for(int i=0;i<Parameter.pathFood*xdim*ydim;i++){
				
				//choose random path to add food to
				if(Math.round(whichPath.nextDouble())==1){
					
					double xpos = 1000 + wherePath.nextDouble();
					double ypos = 1*xpos - 750;
					
					Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
					sitesAdded++;
					allCells.add(cell);
					markerPoint mp = new markerPoint(context,xpos,ypos);
					pathCells.add(mp);
					
				} else {
					
					double xpos = 500 + wherePath.nextDouble();
					double ypos = 1*xpos + 250;
					
					Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
					sitesAdded++;
					allCells.add(cell);
					markerPoint mp = new markerPoint(context,xpos,ypos);
					pathCells.add(mp);
				}
				
			} 
		} else if (Parameter.addPath==4){
			
			//path length
			double pathLength = ydim/(4*2);
			
			//food per path
			double pathFood = Parameter.pathFood*xdim*ydim/4;
			
			cern.jet.random.Normal err = RandomHelper.createNormal(0, Parameter.pathWidth);
			cern.jet.random.Uniform wherePath = RandomHelper.createUniform(0,pathLength);
			
			//Path 1
			for(int i=0;i<pathFood;i++){
				double xpos = 1000 + wherePath.nextDouble();
				double ypos = 1*xpos - 750;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 2
			for(int i=0;i<pathFood;i++){
				double xpos = 500 + wherePath.nextDouble();
				double ypos = 1*xpos + 250;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 3
			for(int i=0;i<pathFood;i++){
				double xpos = 1500 + wherePath.nextDouble();
				double ypos = 1*xpos - 250;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 4
			for(int i=0;i<pathFood;i++){
				double xpos = 250 + wherePath.nextDouble();
				double ypos = 1*xpos + 1250;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
		} else if (Parameter.addPath==8){
			
			//path length
			double pathLength = ydim/(8*2);
			
			//food per path
			double pathFood = Parameter.pathFood*xdim*ydim/8;
			
			cern.jet.random.Normal err = RandomHelper.createNormal(0, Parameter.pathWidth);
			cern.jet.random.Uniform wherePath = RandomHelper.createUniform(0,pathLength);
			
			//Path 1
			for(int i=0;i<pathFood;i++){
				double xpos = 1000 + wherePath.nextDouble();
				double ypos = 1*xpos - 750;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 2
			for(int i=0;i<pathFood;i++){
				double xpos = 500 + wherePath.nextDouble();
				double ypos = 1*xpos + 250;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 3
			for(int i=0;i<pathFood;i++){
				double xpos = 1500 + wherePath.nextDouble();
				double ypos = 1*xpos - 250;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 4
			for(int i=0;i<pathFood;i++){
				double xpos = 250 + wherePath.nextDouble();
				double ypos = 1*xpos + 1250;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 5
			for(int i=0;i<pathFood;i++){
				double xpos = 250 + wherePath.nextDouble();
				double ypos = 1*xpos - 250;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 6
			for(int i=0;i<pathFood;i++){
				double xpos = 1000 + wherePath.nextDouble();
				double ypos = 1*xpos + 750;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 7
			for(int i=0;i<pathFood;i++){
				double xpos = 1750 + wherePath.nextDouble();
				double ypos = 1*xpos - 1250;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
			
			//Path 8
			for(int i=0;i<pathFood;i++){
				double xpos = 1000 + wherePath.nextDouble();
				double ypos = 1*xpos + 0;
				
				Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
				sitesAdded++;
				allCells.add(cell);
				markerPoint mp = new markerPoint(context,xpos,ypos);
				pathCells.add(mp);
			}
		} else if (Parameter.addPath==10){
			
			//path length
			double pathLength = ydim/(1*2);
			
			cern.jet.random.Normal err = RandomHelper.createNormal(0, Parameter.pathWidth);
			cern.jet.random.Uniform wherePath = RandomHelper.createUniform(0,pathLength);
			for(int i=0;i<Parameter.pathFood*xdim*ydim;i++){
				
					double xpos = 1000 + wherePath.nextDouble();
					double ypos = 1*xpos-750;
					
					Cell cell = new Cell(context,xpos+ err.nextDouble(),ypos+ err.nextDouble(),beta.nextDouble(),count++);
					sitesAdded++;
					allCells.add(cell);
					markerPoint mp = new markerPoint(context,xpos,ypos);
					pathCells.add(mp);
				
			} 
		}
		

		/************************************
		 * 							        *
		 * Adding hosts to the landscape	*
		 * 							        *
		 * *********************************/

		//keep track of groups being added
		System.out.println("adding agents to the environment");
		ArrayList<Primate> group = new ArrayList<Primate>();

		//center of group (fixed)
		double xCenter =xdim/2;
		double yCenter =0;


		//select the number of primates in this group (fixed)
		int groupSize = Parameter.groupSize;

		boolean isMale = true;

		//add individuals
		for (int j = 0; j < groupSize; j++){

			//add individual
			Coordinate coord=SimUtils.generateCoordAround(xCenter,yCenter);

			Baboon rc = new Baboon(primatesAdded++,coord,groupSize,isMale);
			isMale=false;
			context.add(rc);
			primatesAll.add(rc);
			orderedP.add(rc);
			group.add(rc);
			Point geom = fac.createPoint(coord);
			geog.move(rc, geom);
			rc.myPatch = null;
		}

		//Add groupMates list (for simulation control)
		for(Baboon p:primatesAll){
			p.setPrimateList(new ArrayList(primatesAll));
		}

		//setup influence structure
		if(Parameter.influenceType==1){
			NetworkUtils.leaderNet(orderedP);
		} else if (Parameter.influenceType==2){
			NetworkUtils.corePeriphery(orderedP,Parameter.corePer);
		} else {
			NetworkUtils.randomNet(orderedP);
		}

		for(Primate p:this.getAllPrimateAgents()){
			System.out.println("I'm "+p.id+ " following " + p.followMate.id);
		}


		/************************************
		 * 							        *
		 * create the scheduling			*
		 * 							        *
		 * *********************************/

		// Ordering processes
		// (1) get inputs for agents						- threaded
		// (2) behavioural responses of agents (movement) 	- threaded
		// (3) environment growback							- threaded
		// (4) add/remove cells to modify					- not threaded

		//executor takes care of the processing of the schedule
		Executor executor = new Executor();
		createSchedule(executor);

		/************************************
		 * 							        *
		 * Adding watcher agent				*
		 * 							        *
		 * *********************************/

		//watcher agent records all output
		Watcher w = new Watcher(executor);
		context.add(w);

		return context;

	}

	private void createSchedule(Executor executor){

		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();

		if(Parameter.groupProcess==true){ //subgroups of agents passed through each separate process, then next subgroup processed

			ScheduleParameters agentStepParamsPrimate = ScheduleParameters.createRepeating(1, 1, 6); //start, interval, priority (high number = higher priority)
			schedule.schedule(agentStepParamsPrimate,executor,"processPrimates");

		} else {  //all agents passed through separate process, one process at a time

			ScheduleParameters agentStepParamsPrimate = ScheduleParameters.createRepeating(1, 1, 6); //start, interval, priority (high number = higher priority)
			schedule.schedule(agentStepParamsPrimate,executor,"getInputs");

			ScheduleParameters agentStepParamsPrimateBehaviour = ScheduleParameters.createRepeating(1, 1, 5); //start, interval, priority (high number = higher priority)
			schedule.schedule(agentStepParamsPrimateBehaviour,executor,"behaviour");

			ScheduleParameters agentStepParamsPrimateEnergy = ScheduleParameters.createRepeating(1, 1, 4); //start, interval, priority (high number = higher priority)
			schedule.schedule(agentStepParamsPrimateEnergy,executor,"energyUpdate");

		}

		ScheduleParameters agentStepParams = ScheduleParameters.createRepeating(1, 1, 2);
		schedule.schedule(agentStepParams,executor,"envUpdate");
	}

	private static void resetP(){

		String csvFile = Parameter.parameters_csv;
		String line = "";
		String cvsSplitBy = ",";
		String[] params_new=null;


		try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {

			while ((line = br.readLine()) != null) {

				// use comma as separator
				params_new = line.split(cvsSplitBy);
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		Parameter.addPath = Integer.parseInt(params_new[0]);
		Parameter.groupSize = Integer.parseInt(params_new[1]);
		Parameter.influenceType = Integer.parseInt(params_new[2]);
		Parameter.corePer = Double.parseDouble(params_new[3]);
		Parameter.pathWidth = Integer.parseInt(params_new[4]);
		Parameter.pathFood = Double.parseDouble(params_new[5]);
		Parameter.stepsPerDay = Double.parseDouble(params_new[6]);


	}



	//used to update only the cells which have been modified, or are growing back
	@ScheduledMethod(start=1, interval = 1,priority=1)
	public static void removeCellsUpdated(){
		for(Cell c : removeCellsToUpdate){
			cellsToUpdate.remove(c);
		}
		removeCellsToUpdate.clear();
	}

	public static void stopSim(Exception ex, Class<?> clazz) {
		ISchedule sched = RunEnvironment.getInstance().getCurrentSchedule();
		sched.setFinishing(true);
		sched.executeEndActions();
	}

	public static Iterable<Cell> getAllCellAgents() {
		return allCells;
	}

	public synchronized static Iterable<Baboon> getAllPrimateAgents(){
		return orderedP;
	}

	public static Context<Cell> getContext() {
		return ModelSetup.mainContext;
	}
	public static Geography<Primate> getGeog(){
		return geog;
	}
	public synchronized static void addToCellUpdateList(Cell c){  //need to read up on the synchronized implimentation
		if(cellsToUpdate.contains(c)==false)cellsToUpdate.add(c);
	}
	public synchronized static void removeCellToUpdate(Cell c){
		removeCellsToUpdate.add(c);
	}
	public static ArrayList<Cell> getCellsToUpdate(){
		return cellsToUpdate;
	}
	public static <T> Geometry getAgentGeometry(T agent) {
		return getGeog().getGeometry(agent);
	}
	public static GeodeticCalculator getGC(){
		return gc;
	}

}
