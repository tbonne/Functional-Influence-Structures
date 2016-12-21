package LHP;

import repast.simphony.visualizationOGL2D.DefaultStyleOGL2D;
import saf.v3d.scene.VSpatial;

public class BaboonStyle2D extends DefaultStyleOGL2D{
	
	/**
	   * @return a circle of radius 4.
	   */
	  public VSpatial getVSpatial(Object agent, VSpatial spatial) {
	    if (spatial == null) {
	      spatial = shapeFactory.createRectangle(2, 5);
	    }
	    return spatial;
	  }
}
