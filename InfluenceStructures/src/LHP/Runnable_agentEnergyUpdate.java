package LHP;

public class Runnable_agentEnergyUpdate implements Runnable{

	Primate primate;

	Runnable_agentEnergyUpdate(Primate p){
		primate = p;
	}

	@Override
	public void run(){
		Throwable thrown = null;
		try {
			//System.out.println("running energy");
			primate.energyUpdate();
		} catch (Throwable e) {
			thrown = e;
			System.out.println("Problem lies in energy update code" + thrown);
		} finally {
			return;
		}
	}

}
