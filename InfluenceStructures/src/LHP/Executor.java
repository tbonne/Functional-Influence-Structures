package LHP;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;


/*
 * 
 * This class controls the multithreading processes, breaking the model up into 5 seperate multithreaded sequences:
 * 			1) getInputs
 * 			2) Behavioural response
 * 			3) Energy update
 * 			4) Memory update
 * 			5) Environment regrowth
 */

public class Executor {

	private static final int pThreads = Parameter.numbOfThreads;
	private static ExecutorService executor;

	public Executor(){
		executor = Executors.newFixedThreadPool(pThreads);
	}

	/*******************************************************************************************************************************/
	
	//n-asynchronous scheduling
	
	public static void processPrimates(){

		//get all primates into an iterator 
		Iterator<Baboon> primateList = ModelSetup.getAllPrimateAgents().iterator();
		
		//System.out.println("started processing primates");
		
		//for all primates in order
		while (primateList.hasNext()){

			//get subsample of primates
			ArrayList<Primate> pToProcess = new ArrayList<Primate>();
			for (int i=0; i<Parameter.numbOfThreads;i++){
				try{
					pToProcess.add(primateList.next());
					//System.out.println("another primate added to group list");
				} catch (java.util.NoSuchElementException e){

				}
			}

			//process subsample of the primate population

			//inputs
			Collection<Callable<Void>> tasks_inputs = new ArrayList<Callable<Void>>();
			for (Primate p:pToProcess){
				Runnable worker = new Runnable_agentInputs(p);
				tasks_inputs.add(Executors.callable(worker,(Void)null));
			}

			try {
				for (Future<?> f : executor.invokeAll(tasks_inputs)) { //invokeAll() blocks until ALL tasks submitted to executor complete
					f.get(); 
					//System.out.println("inputs " + pThreads + Parameter.groupProcess);
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}catch (NullPointerException e){
				e.printStackTrace();
			}

			//behaviour
			Collection<Callable<Void>> tasks_behaviour = new ArrayList<Callable<Void>>();
			for (Primate p:pToProcess){
				Runnable worker = new Runnable_agentBehaviour(p);
				tasks_behaviour.add(Executors.callable(worker,(Void)null));
			}

			try {
				for (Future<?> f : executor.invokeAll(tasks_behaviour)) { //invokeAll() blocks until ALL tasks submitted to executor complete
					f.get();
					//System.out.println("behaviour " + pThreads + Parameter.groupProcess);
				}
				
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}catch (NullPointerException e){
				e.printStackTrace();
			}
			
			//energy uptade
			Collection<Callable<Void>> tasks_energy = new ArrayList<Callable<Void>>();
			for (Primate p:pToProcess){
				Runnable worker = new Runnable_agentEnergyUpdate(p);
				tasks_energy.add(Executors.callable(worker,(Void)null));
			}

			try {
				for (Future<?> f : executor.invokeAll(tasks_energy)) { //invokeAll() blocks until ALL tasks submitted to executor complete
					f.get(); 
					//System.out.println("energy");
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}catch (NullPointerException e){
				e.printStackTrace();
			}
			
		} //loop: process new subsample

	}

	/*******************************************************************************************************************************/
	
	//Synchronous scheduling
	
	public static void getInputs(){
		//System.out.println("Starting input threads");

		Collection<Callable<Void>> tasks = new ArrayList<Callable<Void>>();
		Iterator<Baboon> primateList = ModelSetup.getAllPrimateAgents().iterator();
		while (primateList.hasNext()){
			Runnable worker = new Runnable_agentInputs(primateList.next());
			tasks.add(Executors.callable(worker,(Void)null));
		}

		// Wait until all threads are finish
		try {
			for (Future<?> f : executor.invokeAll(tasks)) { //invokeAll() blocks until ALL tasks submitted to executor complete
				f.get(); 
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		} catch (NullPointerException e){
			e.printStackTrace();
			System.out.println("null in inputs");
		}

		//System.out.println("Finished all input threads");

	}

	public static void behaviour(){
		Collection<Callable<Void>> tasks = new ArrayList<Callable<Void>>();
		Iterator<Baboon> primateList = ModelSetup.getAllPrimateAgents().iterator();
		while (primateList.hasNext()){
			Runnable worker = new Runnable_agentBehaviour(primateList.next());
			tasks.add(Executors.callable(worker,(Void)null));
		}

		// Wait until all threads are finish
		try {
			for (Future<?> f : executor.invokeAll(tasks)) { //invokeAll() blocks until ALL tasks submitted to executor complete
				f.get(); 
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		} catch (NullPointerException e){
			e.printStackTrace();
			System.out.println("null in behaviour");
		}

		//System.out.println("Finished all behaviour threads");
	}


	public static void energyUpdate(){
		Collection<Callable<Void>> tasks = new ArrayList<Callable<Void>>();
		Iterator<Baboon> primateList = ModelSetup.getAllPrimateAgents().iterator();
		while (primateList.hasNext()){
			Runnable worker = new Runnable_agentEnergyUpdate(primateList.next());
			tasks.add(Executors.callable(worker,(Void)null));
		}

		// Wait until all threads are finish
		try {
			for (Future<?> f : executor.invokeAll(tasks)) { //invokeAll() blocks until ALL tasks submitted to executor complete
				f.get(); 
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		} catch (NullPointerException e){
			e.printStackTrace();
			System.out.println("null in energy updates");
		}

		//System.out.println("Finished all energy threads");
	}

	public static void envUpdate(){
		Collection<Callable<Void>> tasks = new ArrayList<Callable<Void>>();

		Iterable<Cell> cells = ModelSetup.getCellsToUpdate();
		for (Cell c: cells){
			Runnable worker = new Runnable_envUpdate(c);
			tasks.add(Executors.callable(worker,(Void)null));
		}

		// Wait until all threads are finish
		try {
			for (Future<?> f : executor.invokeAll(tasks)) { //invokeAll() blocks until ALL tasks submitted to executor complete
				f.get(); 
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}catch (NullPointerException e){
			e.printStackTrace();
		}
		//System.out.println("Finished all cell threads");
	}

	public void shutdown(){
		executor.shutdown();
	}
}
