// STAR-CCM+ macro: Arearule.java
//Method calculates the cross sectional area at each given x location along the aircraft
//Can be used for any mach number
//The model needs to be set up so that the freestream is a cylinder around the model that is wider and longer than the model
//Outputs in the log file
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.vis.*;
import star.base.report.*;

public class AreaRule extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

	//Mach number to evaluate area rule for
	double m = 1.8;

	
	int i = 0;
	
	//Set up constants for mach angle
	Double totalArea = 0.0;
	double mu = Math.asin(1.0/m);
	double sin = Math.sin(mu);
	double cos = Math.cos(mu);
	simulation_0.println(mu +"," +sin + "," +cos);
	
	//Iterate along x axis of aircraft
	//for (Double x = -10.0; x < 50.0; x += 1.0) {
	for (Double x = 0.; x < 60.; x += 1.0) { 

 		//Grab Plane
		PlaneSection planeSection_0 = 
			((PlaneSection) simulation_0.getPartManager().getObject("plane section"));
   
		Coordinate coordinate_0 = 
			planeSection_0.getOriginCoordinate();

		Units units_0 = 
			((Units) simulation_0.getUnitsManager().getObject("m"));

		//Set origin for plane
		coordinate_0.setCoordinate(units_0, units_0, units_0, new DoubleVector(new double[] {x, 0.0, 0.0}));
		
		//Counting variables for cross sectional area
		double runningTotal = 0.0;
		double total = 0.0;
		//Denotes the number of planes at each x location- Use more for higher accuracy
		double iter = 4.0;
		
		//Iterate through planes at x location
		for (double t = 0.0; t < 2.0 * Math.PI; t += Math.PI / iter ) {

			//Set plane orientation
    			Coordinate coordinate_1 = 
      				planeSection_0.getOrientationCoordinate();


    			coordinate_1.setCoordinate(units_0, units_0, units_0, new DoubleVector(new double[] {sin, Math.cos(t) * cos, Math.sin(t) * cos}));

    
			//Calculate frontal area
    			FrontalAreaReport frontalAreaReport_0 = 
      				((FrontalAreaReport) simulation_0.getReportManager().getReport("Frontal Area 1"));

    			//Add up areas
			Double fieldValue = frontalAreaReport_0.getReportMonitorValue();
			runningTotal+=fieldValue;
			
		}
		
		//Average areas
		runningTotal = runningTotal / iter;
		
		//Find the initial offset area
		if (i == 0) {
			totalArea = runningTotal;	
		}

		//Scale area due to plane angle
		double value = (totalArea - runningTotal) * sin;
		//simulation_0.println(x + "," + totalArea - fieldValue.doubleValue());
		//Print area
		simulation_0.println(x + "," + value);
		i++;
		
	}


  }
}
