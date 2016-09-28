package slather.g7;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;

import slather.sim.Cell;
import slather.sim.GridObject;
import slather.sim.Move;
import slather.sim.Pherome;
import slather.sim.Point;

public class Player implements slather.sim.Player {

	private Random gen;

	@Override
	public void init(double d, int t) {
		gen = new Random();
	}

	@Override
	public Move play(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {

		// Convert the byte memory to binary string for easier usage
		Memory m = new Memory(memory);

		/*
		 * If cell can reproduce, do it. Otherwise, move based on clustering
		 * strategy
		 */

		// Reproduction
		if (player_cell.getDiameter() >= 2) {
			byte childMemory = generateChildMemories(memory);
			return new Move(true, childMemory, oppositeChild(childMemory).getByte());
		}

		Point moveDirection = getCumulativeDirection(player_cell, nearby_cells, nearby_pheromes);

		// Last bit denotes direction of movement
		// If it is 1, move towards our own cells
		if (m.opposite == 1) {
			moveDirection = new Point(-moveDirection.x, -moveDirection.y);
		}

		Point finalMoveDirection = addRandomnessAndGenerateFinalDirection(moveDirection);

		byte nextMoveMemory = generateNextMoveMemory(memory);

		return new Move(finalMoveDirection, nextMoveMemory);
	}

	/* Scaffolding zone: Fill in the functions here */
	public byte generateChildMemories(byte memory) {
		// Child will go only 2 steps away for now
		Memory child = new Memory(3, 2, 0);

		return child.getByte();
	}

	public Point getCumulativeDirection(Cell myCell, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
		System.out.println("Join forces from cells:");
		Point fromCells = joinForcesFromCells(myCell, nearby_cells, 1.0);
		System.out.println("Join forces from pheromes:");
		Point fromPheromes = joinForcesFromPheromes(myCell, nearby_pheromes, 0.5);
		Point normalized = normalizeDistance(fromCells, fromPheromes);

		return normalized;
	}

	private Point addRandomnessAndGenerateFinalDirection(Point direction) {
		double desiredMean = 0.0;
		double desiredStandardDeviation = 0.1;// Set a very small deviation so
												// that the whole force thing
												// still makes sense
		double offX = direction.x + gen.nextGaussian() * desiredStandardDeviation + desiredMean;
		double offY = direction.y + gen.nextGaussian() * desiredStandardDeviation + desiredMean;

		Point newDir = normalizeDistance(new Point(offX, offY));
		System.out.println("The randomized direction has X: " + newDir.x + " Y: " + newDir.y);
		return newDir;
	}

	/* Helper functions */
	private byte generateNextMoveMemory(byte memory) {
		Memory memoryObject = new Memory(memory);
		if (memoryObject.opposite == 1) {
			
			if (memoryObject.moveDirectionCountdown == 0) {
				memoryObject.opposite = 0;
				memoryObject.offsetCountDown = 3;
				
				// Setting countdown based on cell life
				if (memoryObject.moveDirectionCountSetter >= 2) {
					memoryObject.moveDirectionCountSetter = 3;
					memoryObject.moveDirectionCountdown = 7;
				} else {
					(memoryObject.moveDirectionCountSetter)++;
					memoryObject.moveDirectionCountdown = 2 * memoryObject.moveDirectionCountSetter;
				}
				
			} else {
				memoryObject.moveDirectionCountdown--;
			}
			
		} else {
			
			if (memoryObject.moveDirectionCountdown == 0) {
				memoryObject.opposite = 1;
				memoryObject.offsetCountDown = 0;

				// Setting countdown based on cell life
				if (memoryObject.moveDirectionCountSetter >= 2) {
					memoryObject.moveDirectionCountSetter = 3;
					memoryObject.moveDirectionCountdown = 7;
				} else {
					(memoryObject.moveDirectionCountSetter)++;
					memoryObject.moveDirectionCountdown = 2 * memoryObject.moveDirectionCountSetter;
				}
				
			} else {
				if (memoryObject.offsetCountDown == 0) {
					memoryObject.moveDirectionCountdown--;
				} else {
					memoryObject.offsetCountDown--;
				}
			}
			
		}

		return memoryObject.getByte();
	}

	/*
	 * Based on the assumption that we chase opponent cells and get away from
	 * friendly ones.
	 */
	public static Point joinForcesFromCells(Cell myCell, Set<Cell> cells, double weight) {
		Point myPos = myCell.getPosition();

		double offX = 0.0;
		double offY = 0.0;
		for (Cell c : cells) {
			Point cPos = c.getPosition();
			/* The smaller the distance, the larger the force */
			int chaseIt;
			if (c.player == myCell.player)
				chaseIt = -1;
			else
				chaseIt = 1;

			double diffX = cPos.x - myPos.x;
			double diffY = cPos.y - myPos.y;
			// The force is proportional to the mass of the cell, which depends
			// on the square of diameter in a 2d environment
			offX += chaseIt * diffX * weight * (c.getDiameter() * c.getDiameter());
			offY += chaseIt * diffY * weight * (c.getDiameter() * c.getDiameter());
		}

		System.out.println("The merged force from cells is X: " + offX + " Y: " + offY);

		return new Point(offX, offY);
	}

	public static Point joinForcesFromPheromes(Cell myCell, Set<Pherome> pheromes, double weight) {
		Point myPos = myCell.getPosition();
		double offX = 0.0;
		double offY = 0.0;
		for (Pherome p : pheromes) {
			Point pPos = p.getPosition();
			/* The smaller the distance, the larger the force */
			int chaseIt;
			if (p.player == myCell.player)
				chaseIt = -1;
			else
				chaseIt = 1;

			double diffX = pPos.x - myPos.x;
			double diffY = pPos.y - myPos.y;
			// The force is proportional to the mass of the cell, which depends
			// on the square of diameter in a 2d environment
			offX += chaseIt * diffX * weight;
			offY += chaseIt * diffY * weight;
		}

		System.out.println("The merged force from pheromes is X: " + offX + " Y: " + offY);
		return new Point(offX, offY);
	}

	public static Point normalizeDistance(Point... points) {
		double offX = 0.0;
		double offY = 0.0;
		for (Point p : points) {
			offX += p.x;
			offY += p.y;
		}
		// Normalize the force to 1mm length
		double hypotenuse = Math.sqrt(offX * offX + offY * offY);

		Point toReturn;
		if (hypotenuse == 0.0) {
			toReturn = new Point(0.0, 0.0);
		} else {
			offX /= hypotenuse;
			offY /= hypotenuse;
			toReturn = new Point(offX, offY);
		}

		System.out.println("The normalized direction is X: " + toReturn.x + " Y: " + toReturn.y);
		return toReturn;
	}

	public static void reportMemory(Memory m) {
		System.out.println("Memory has offsetCountDown: " + m.offsetCountDown);
		System.out.println("Memory directionCounter:" + m.moveDirectionCountdown);
		System.out.println("Whether to go opposite:" + m.opposite);

	}

	private boolean collides(Cell player_cell, Point vector, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
		Iterator<Cell> cell_it = nearby_cells.iterator();
		Point destination = player_cell.getPosition().move(vector);
		while (cell_it.hasNext()) {
			Cell other = cell_it.next();
			if (destination.distance(other.getPosition()) < 0.5 * player_cell.getDiameter() + 0.5 * other.getDiameter()
					+ 0.00011)
				return true;
		}
		Iterator<Pherome> pherome_it = nearby_pheromes.iterator();
		while (pherome_it.hasNext()) {
			Pherome other = pherome_it.next();
			if (other.player != player_cell.player
					&& destination.distance(other.getPosition()) < 0.5 * player_cell.getDiameter() + 0.0001)
				return true;
		}
		return false;
	}

	public static Memory oppositeChild(Memory m) {
		return new Memory(m.offsetCountDown, m.moveDirectionCountdown, 1 - m.opposite);
	}

	public static Memory oppositeChild(byte b) {
		Memory m = new Memory(b);
		return new Memory(m.offsetCountDown, m.moveDirectionCountdown, 1 - m.opposite);
	}

}
