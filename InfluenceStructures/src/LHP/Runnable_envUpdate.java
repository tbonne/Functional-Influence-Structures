package LHP;

public class Runnable_envUpdate implements Runnable{
	
	Cell cell;
	
	Runnable_envUpdate(Cell c){
		cell = c;
	}
	
	@Override
	public void run(){
		Throwable thrown = null;
	    try {
		cell.stepThreaded();
	    } catch (Throwable e) {
	        thrown = e;
	        System.out.println("Problem lies in cell update: "+ thrown);
	    } finally {
	        //threadExited(this, thrown);
	    	return;
	    }
	}

}
