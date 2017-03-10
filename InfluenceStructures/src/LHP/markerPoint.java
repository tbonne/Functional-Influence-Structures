package LHP;

import java.util.ArrayList;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.space.gis.Geography;

public class markerPoint {

	//geometry
	Context context;
	Geography geog;
	Geometry geom;
	Coordinate coord; //centroid

	public markerPoint(Context con,double x,double y){

		//set initial variables
		con.add(this);

		//place the cell on a gis landscape
		geog = (Geography)con.getProjection("geog");
		geom = getCircle(x,y);
		//geom = getHexShape(x,y);

		geog.move(this, geom);
		this.setCoord(geom.getCentroid().getCoordinate());

	}

	private Geometry getCircle(double x, double y){

		GeometryFactory fac = new GeometryFactory();
		Geometry shape = null;
		Coordinate centroid = new Coordinate(x,y);

		Point point = fac.createPoint(centroid);
		shape = point.buffer(Parameter.foodBuffer);

		return shape;
	}

	private void setCoord(Coordinate c){
		coord = c;
	}
	public Coordinate getCoord(){
		return this.coord;
	}
}
