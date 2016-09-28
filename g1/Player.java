package slather.g1;

import slather.sim.Cell;
import slather.sim.Move;
import slather.sim.Pherome;
import slather.sim.Point;

import java.util.*;

public class Player implements slather.sim.Player {

	private int playerType;
	private double distanceVisible;
	private static final double MAXIMUM_MOVE = 1.0;
	private static final double THRESHOLD_DISTANCE = 2.0;
	private Random gen;

	public void init(double d, int t) {
		playerType = t;
		distanceVisible = d;
		gen = new Random();
	}

	public Move play(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
		if (player_cell.getDiameter() >= 2) // reproduce whenever possible
			return new Move(true, (byte)-1, (byte)-1);

		List<Vector> vectors = new ArrayList<>();
		//Get list of vectors of nearby neighbors
		Point myPosition = player_cell.getPosition();
		for(Cell c : nearby_cells) {
			Point otherPosition = c.getPosition();

			//Check distance between cell and neighboring cell. If greater than theshold, add vector
//			if(myPosition.distance(otherPosition) < THRESHOLD_DISTANCE)
				vectors.add(new Vector(myPosition, otherPosition));
		}
		//Also get list of nearby pheremones
//		for(Pherome p : nearby_pheromes) {
//			Point otherPosition = p.getPosition();
//
//			//Check distance between cell and neighboring cell. If greater than theshold, add vector
////			if(myPosition.distance(otherPosition) < THRESHOLD_DISTANCE)
//			if(playerType != p.player)
//				vectors.add(new Vector(myPosition, otherPosition));
//		}

		//Add all vectors together
		Vector finalVector = new Vector(0, 0);
		for(Vector v : vectors) {
			finalVector = finalVector.add(v);
		}

		//Get inverse of direction
		finalVector = finalVector.invert();

		//Get final destination point
		Point finalPoint = finalVector.add(myPosition);

		//Make sure point falls within MAXIMUM_MOVE
		double distance = finalPoint.distance(myPosition);
		if(distance > MAXIMUM_MOVE) {
			finalPoint = finalVector.multiply(MAXIMUM_MOVE/distance).add(myPosition);
		}

//		System.out.println("Old " + myPosition);
//		System.out.println("Final Point:" + finalPoint);
//		System.out.println("Old Distance " + distance);
//		System.out.println("New Distance " + myPosition.distance(finalPoint));
//		System.out.println();

		// if all tries fail, just chill in place
		return new Move((new Vector(myPosition, finalPoint)).toPoint(), (byte)0);
	}

}
