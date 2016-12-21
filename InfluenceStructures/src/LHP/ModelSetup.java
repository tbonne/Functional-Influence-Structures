package LHP;

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
import tools.SimUtils;

//This file builds the model: creating the environment and populating it with primates
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
		mainContext = context; //static link to context
		timeRecord = System.currentTimeMillis();
		timeRecord_start = System.currentTimeMillis();

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


		//add Resources to the environment (high density path)
/*		Uniform xDist2 = RandomHelper.createUniform(xdim/2, xdim);
		cern.jet.random.Normal err = RandomHelper.createNormal(0, 10);
		for(int i=0;i<Parameter.foodDensity*xdim*ydim;i++){
			double xSample = xDist2.nextDouble();
			double ySample = ydim+1;
			while (ySample > ydim || ySample < 0){
				//ySample = -1*Math.pow((xSample-xdim/2)/2,2)+ydim/2+err.nextDouble();
				ySample = Math.pow( (xSample-xdim/2.0), 0.25) ;
				ySample = ySample / Math.pow(xdim/2,0.25)          ;
				ySample = ySample * ydim ;
			}

			Cell cell = new Cell(context,xSample+ err.nextDouble(),ySample+ err.nextDouble(),beta.nextDouble(),count++);
			sitesAdded++;
			allCells.add(cell);
		}
*/
		//add Resources to the environment (high density path)
		Uniform yDist2 = RandomHelper.createUniform(0, ydim);
		cern.jet.random.Normal err = RandomHelper.createNormal(0, 10);
		for(int i=0;i<Parameter.foodDensity*xdim*ydim/4;i++){
			double ySample = yDist2.nextDouble();
			double xSample = xdim+1;
			while (xSample > xdim || xSample < 0){
				//ySample = -1*Math.pow((xSample-xdim/2)/2,2)+ydim/2+err.nextDouble();
				xSample = Math.pow( ySample, 0.25) ;
				xSample = xSample / Math.pow(ydim,0.25)          ;
				xSample = xSample * xdim ;
			}

			Cell cell = new Cell(context,xSample+ err.nextDouble(),ySample+ err.nextDouble(),beta.nextDouble(),count++);
			sitesAdded++;
			allCells.add(cell);
		}


		/*
		//adding hexagon grid cells
		Beta beta = new Beta(Parameter.envHomogen,Parameter.envHomogen, Parameter.mt);
		double xcoord=0,ycoord=0;
		int offset=0,count=0;

		for (int i = 0; i < ydim; ++i) {
			for (int j = 0; j < xdim; ++j) {

				//double food = 0;
				double food = beta.nextDouble();
				Cell cell=null;

				if(offset==0){
					xcoord=xcoord+(Parameter.cellSize/2.0)*Math.cos(2*Math.PI/3.0);
					offset=1;
					cell = new Cell(context,xcoord,ycoord,food,count++);
					xcoord=xcoord-(Parameter.cellSize/2.0);
				}else if (offset==1){
					xcoord=xcoord-(Parameter.cellSize/2.0)*Math.cos(2*Math.PI/3.0);
					ycoord=ycoord+(Parameter.cellSize/2.0)*Math.sin(2*Math.PI/3.0);
					offset=0;
					cell = new Cell(context,xcoord,ycoord,food,count++);
					ycoord=ycoord-(Parameter.cellSize/2.0)*Math.sin(2*Math.PI/3.0);
				}

				resAdded=resAdded+food;
				allCells.add(cell);

				//shift ycoord by cell size value
				xcoord=xcoord+Parameter.cellSize;
			}

			//Set ycoord back to the start and shift xcoord up by cell size value
			xcoord=0;
			ycoord = ycoord-(Parameter.cellSize)*Math.sin(2*Math.PI/3.0);
		}
		 */

		//to simplify the model all cells record the visible neighbours (within visible range for a primate)
		for (Cell c: allCells){
			c.setVisibleNeigh();
		}

		/************************************
		 * 							        *
		 * Adding hosts to the landscape	*
		 * 							        *
		 * *********************************/

		//keep track of groups being added
		ArrayList<Primate> group = new ArrayList<Primate>();

		//center of group (fixed)
		double xCenter =xdim/2;//*Parameter.cellSize/3;//0;//xDist.nextDouble();//75;//1000;// 75;
		double yCenter =0;//0;//yDist.nextDouble();// -75;//-1000;//-75;
		//double yoffset = 10;
		//double xoffset = 0;

		for(int i = 0 ;i<Parameter.numbOfGroups;i++ ){

			//select the number of primates in this group (fixed)
			int groupSize = Parameter.groupSize;

			boolean isMale = true;

			//add individuals
			for (int j = 0; j < groupSize; j++){

				//add individual
				Coordinate coord=SimUtils.generateCoordAround(xCenter,yCenter);
				//Coordinate coord= new Coordinate(xCenter+RandomHelper.nextDouble()*0.001,yCenter+RandomHelper.nextDouble()*0.001);
				//xCenter = xCenter + xoffset;
				//yCenter = yCenter + yoffset;

				Baboon rc = new Baboon(primatesAdded++,coord,groupSize,isMale,i);
				isMale=false;
				context.add(rc);
				primatesAll.add(rc);
				orderedP.add(rc);
				group.add(rc);
				Point geom = fac.createPoint(coord);
				geog.move(rc, geom);
				rc.myPatch = null;
			}
		}

		//Add groupMates list (for simulation control)
		for(Baboon p:primatesAll){
			p.setPrimateList(new ArrayList(primatesAll));
		}
		
		for(Primate p:orderedP){
			//System.out.println("starting with "+p.id);
			if(p.getId()==0)p.setBaboonFollower(0);
			if(p.getId()==1)p.setBaboonFollower(0);
			if(p.getId()==2)p.setBaboonFollower(0);//3
			if(p.getId()==3)p.setBaboonFollower(0);
			if(p.getId()==4)p.setBaboonFollower(2);
			if(p.getId()==5)p.setBaboonFollower(2);
			if(p.getId()==6)p.setBaboonFollower(2);
			if(p.getId()==7)p.setBaboonFollower(1);
			if(p.getId()==8)p.setBaboonFollower(1);//7
			if(p.getId()==9)p.setBaboonFollower(1);
			if(p.getId()==10)p.setBaboonFollower(3);
			if(p.getId()==11)p.setBaboonFollower(3);
			if(p.getId()==12)p.setBaboonFollower(3);
			//if(p.getId()==13)p.setBaboonFollower(11);
			
			//p.setBaboonFollowerRand(p);
			//p.setBaboonFollowerLeader();
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
