package LHP;

import java.util.ArrayList;
import java.util.Random;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.space.gis.Geography;
import tools.SimUtils;

public class Cell {
	
	//geometry
	Context context;
	Geography geog;
	Geometry geom;
	Coordinate coord; //centroid
	
	//resources
	double resources;
	double maxResources;
	double regrowRate =Parameter.regrowthRate;
	double timeCounter=0;
	
	//memory
	int remembered =0;
	
	//visit counts
	double visit=0;
	int id;
	
	//movement influence
	ArrayList<Cell> visibleSites;
	ArrayList<Double> familiarity;
	
	public Cell(Context con,double x,double y, double r, int i){
		
		//set initial variables
		id=i;
		resources = r;
		maxResources = resources; //maintains initial conditions
		con.add(this);
		
		//place the cell on a gis landscape
		geog = (Geography)con.getProjection("geog");
		geom = getCircle(x,y);
		//geom = getHexShape(x,y);
		
		geog.move(this, geom);
		this.setCoord(geom.getCentroid().getCoordinate());
		timeCounter=0;
		
		visibleSites = new ArrayList<Cell>();
		familiarity = new ArrayList<Double>(Parameter.numbOfGroups*Parameter.groupSize);
		for(int n =0;n<Parameter.numbOfGroups*Parameter.groupSize;n++){
			familiarity.add(0.001);
		}
	}
	
	/*private Geometry getHexShape(double x, double y){

		GeometryFactory fac = new GeometryFactory();
		Geometry shape = null;
		Coordinate[] boundingCoords = new Coordinate[7];
		for(int i=0; i<7; i++) {
			if(i!=6){
				Coordinate c = new Coordinate(x + (Parameter.cellSize/(2.0))*Math.cos(i*2*Math.PI/6), y + (Parameter.cellSize/(2.0))*Math.sin(i*2*Math.PI/6));
				boundingCoords[i]=c;
			} else{
				Coordinate c = new Coordinate(x + (Parameter.cellSize/(2.0))*Math.cos(0*2*Math.PI/6), y + (Parameter.cellSize/(2.0))*Math.sin(0*2*Math.PI/6));
				boundingCoords[i]=c;
			}
		}
		LinearRing ring = fac.createLinearRing(boundingCoords);
		shape = fac.createPolygon(ring, null);

		return shape;
	}*/
	
	private Geometry getCircle(double x, double y){
		
		GeometryFactory fac = new GeometryFactory();
		Geometry shape = null;
		Coordinate centroid = new Coordinate(x,y);
		
		Point point = fac.createPoint(centroid);
		shape = point.buffer(Parameter.foodBuffer);
		
		return shape;
	}
	
	public void setVisibleNeigh(){
		Envelope envelope = new Envelope();
		envelope.init(coord);
		envelope.expandBy(Parameter.visual_range);
		Iterable neigh = geog.getObjectsWithin(envelope, Cell.class);
		while (neigh.iterator().hasNext()){
			Cell neighC = (Cell)neigh.iterator().next();
			if(neighC.getCoord().distance(this.getCoord())<Parameter.visual_range){
				visibleSites.add(neighC);	
			}
			
		}
		envelope.setToNull();
	}
	
	
	//@ScheduledMethod(start=0, interval = 1,priority=2,shuffle=true)
	public void stepThreaded(){
		
		try{
		regrow();
		change();
		}catch (NullPointerException e){
			System.out.println("problem in the regrow method");
		}
		
	}
	
	private void regrow(){
		
		//check to see if resource is at max levels
		if (resources!=maxResources){
			//the resource will grow back until it's maximum level is reached
			if (resources<maxResources-Parameter.regrowthRate){
				resources = (resources+Parameter.regrowthRate);
			} else {
				resources = (maxResources);
			}
			//will degenerate if above max resource level (seasonality)
			if(resources>maxResources+Parameter.regrowthRate){
				resources = (resources-Parameter.regrowthRate);
			}else if (resources>maxResources && resources<maxResources+Parameter.regrowthRate){
				resources = maxResources;
			}
		} else {
			ModelSetup.removeCellToUpdate(this);
		}
	}
	
	private void change(){
		
		for(int i =0;i<familiarity.size();i++){
			double f = familiarity.get(i);
				if(f<0.001){
					f=0.001;
				} else {
					f = Math.max( (f - Parameter.cellChangeRate * (1-f)*f),0.001);
					familiarity.set(i, f);
				}
			}
		}
	
	public synchronized double eatMe(double bite){
		double biteSize =0;
		if (this.getResourceLevel() - bite > 0){
				setResourceLevel( (this.getResourceLevel() - bite));
				biteSize = bite;
		}else{
			biteSize = this.getResourceLevel();
			setResourceLevel(0);
		}
		
		//biteSize = this.getResourceLevel();
		//setResourceLevel(0);
		
		ModelSetup.addToCellUpdateList(this);
		timeCounter=0;
		return biteSize;
	}
	
	/****************************Get and set methods************************************/
	
	public double getResourceLevel(){
		return resources;
	}
	public void setResourceLevel(double r){
		resources = r;
	}
	public double getMaxResourceLevel(){
		return maxResources;
	}
	public void setMaxResourceLevel(double rm){
		maxResources = rm;
	}
	public Coordinate getCoord(){
		return coord;
	}
	private void setCoord(Coordinate c){
		coord = c;
	}
	public int getID(){
		return id;
	}
}
