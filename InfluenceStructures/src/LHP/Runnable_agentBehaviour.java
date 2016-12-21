package LHP;

import java.util.Arrays;

public class Runnable_agentBehaviour implements Runnable {


	Primate primate;

	Runnable_agentBehaviour(Primate p){
		primate = p;
	}

	@Override
	public void run(){
		Throwable thrown = null;
		try {
			//System.out.println("running behaviour");
			primate.behaviouralResponse();
		} catch (Throwable e) {
			thrown = e;
			System.out.println("Problem lies in behaviour code" + thrown);
			//System.out.println(org.apache.commons.lang.exception.ExceptionUtils.getFullStackTrace(e));
		} finally {
			return;
		}
	}


}
