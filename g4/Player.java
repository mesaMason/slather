package slather.g4;

import slather.sim.Cell;
import slather.sim.Point;
import slather.sim.Move;
import slather.sim.Pherome;
import java.util.*;

public class Player implements slather.sim.Player {

	private Random gen;

	public void init(double d, int t) {
		gen = new Random();
	}

	public Move play(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
		if (player_cell.getDiameter() >= 2) // reproduce whenever possible
			return new Move(true, (byte) -1, (byte) -1);
		/* next step
		if (nearby_cells.size() == 1) {
			for (Cell c : nearby_cells) {
				if (c.player == player_cell.player) {
					// System.out.println("1---------------------------------------");
					double newX = player_cell.getPosition().x - c.getPosition().x;
					double newY = player_cell.getPosition().y - c.getPosition().y;
					Point pos = new Point(newX, newY);
					int arg = (int) (Math.atan2(newY, newX) / 3.1415926 * 180) + (newX < 0 ? 180 : 0);
					// System.out.println(arg +
					// "-----------------------------------------");
					Point vector = extractVectorFromAngle(arg);

					if (!collides(player_cell, vector, nearby_cells, nearby_pheromes)) {

						return new Move(vector, (byte) arg);
					}

				}
			}
		}
		*/

		if (nearby_cells.size() == 0) {
			Point vector = extractVectorFromAngle((int) memory);
			// check for collisions
			if (!collides(player_cell, vector, nearby_cells, nearby_pheromes))
				return new Move(vector, memory);
		}

		double[] direction = new double[4];

		for (Cell nearby_cell : nearby_cells) {
			double m_x = nearby_cell.getPosition().x - player_cell.getPosition().x;
			double m_y = nearby_cell.getPosition().y - player_cell.getPosition().y;
			if (m_x >= 0 && m_y >= 0) {
				direction[0] += trans1(nearby_cell.distance(player_cell));
			}
			if (m_x >= 0 && m_y < 0) {
				direction[3] += trans1(nearby_cell.distance(player_cell));
			}
			if (m_x < 0 && m_y >= 0) {
				direction[1] += trans1(nearby_cell.distance(player_cell));
			}
			if (m_x < 0 && m_y < 0) {
				direction[2] += trans1(nearby_cell.distance(player_cell));
			}
		}

		double sum = 0;

		for (int index = 0; index < direction.length; ++index) {
			sum += trans2(direction[index]);
		}

		// We give each direction a probability direction[i]/sum
		sum *= Math.random();
		int index = 0;
		for (; index < direction.length; ++index) {
			sum -= trans2(direction[index]);
			if (sum <= 0) break;
		}

		// int arg = gen.nextInt(30) + min_index * 90 + 30;
		int arg = index * 90 + gen.nextInt(90);

		Point vector = extractVectorFromAngle(arg);
		if (!collides(player_cell, vector, nearby_cells, nearby_pheromes)) {
			return new Move(vector, (byte) arg);
		}
		for (int i = 0; i < 4; i++) {
			int arg2 = gen.nextInt(360) + 1;
			Point vector2 = extractVectorFromAngle(arg2);
			if (!collides(player_cell, vector2, nearby_cells, nearby_pheromes))
				return new Move(vector2, (byte) arg2);
		}
		return new Move(new Point(0, 0), (byte) 0);	
	}

	private double trans1(double a) {
		return 1.0/(a+0.5);
	}

	private double trans2(double a) {
		return 1/Math.pow((0.01+a),2);
	}

	// check if moving player_cell by vector collides with any nearby cell or
	// hostile pherome
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

	// convert an angle (in 2-deg increments) to a vector with magnitude
	// Cell.move_dist (max allowed movement distance)
	private Point extractVectorFromAngle(int arg) {
		double theta = Math.toRadians((double) arg);
		double dx = Cell.move_dist * Math.cos(theta);
		double dy = Cell.move_dist * Math.sin(theta);
		return new Point(dx, dy);
	}

}
