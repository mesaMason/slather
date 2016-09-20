package slather.g8;

import slather.sim.Cell;
import slather.sim.Point;
import slather.sim.Move;
import slather.sim.Pherome;
import java.util.*;


public class Player implements slather.sim.Player {

    private Random gen;
    private double d;
    private int t;

    private static final int NUMDIRECTIONS = 4;   // CONSTANTS - number of directions
    private static final int DURATION_CAP = 4;   // CONSTANTS - cap on traveling
    private static final int FRIENDLY_AVOID_DIST = 5; // CONSTANTS - for random walker: turn away from friendlies that are in this radius
    
    public void init(double d, int t) {
	gen = new Random();
	this.d = d;
	this.t = t;
    }

    /*
      Memory byte setup:
      1 bit strategy (currently either square walker or random walker)

      Memory byte for square walker:
      1 bit strategy
      3 bits empty for now
      2 bits duration
      2 bits direction
      
      Memory byte for random walker:
      1 bit strategy
      7 bits previous direction

      Directions:
      0 = North
      1 = East
      2 = South
      3 = West
     */
    public Move play(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {

        Point currentPosition  = player_cell.getPosition();
        int direction = getDirection(memory);
        int duration = getDuration(memory);
        int maxDuration = getMaxDuration(t, NUMDIRECTIONS);
        int strategy = getStrategy(memory);

	if (player_cell.getDiameter() >= 2) // reproduce whenever possible
	    return new Move(true, (byte) 0b10000000, (byte) 0);
        
        if (strategy == 0) {
            /* Square walker strategy
             */
            duration++;
            if (duration == maxDuration) {
                direction++;
                direction = direction % NUMDIRECTIONS;
                memory = writeDirection(direction, memory);
            }
            memory = writeDuration(duration, memory, maxDuration);

            Point destination = getNewDest(direction);

            /*
            String s = String.format("%8s", Integer.toBinaryString(memory & 0xFF)).replace(' ','0');
            String p = String.format("current x is %f and current y is %f", currentPosition.x, currentPosition.y);
            String des = String.format("dest x is %f and dest y is %f", destination.x, destination.y);
            System.out.println(s);
            System.out.println(p);
            System.out.println(des);
            */
            
            // if collides, try 20 different random directions
            if (collides(player_cell, destination, nearby_cells, nearby_pheromes)) {
                for (int i = 0; i < 20; i++) {
                    int arg = gen.nextInt(120);
                    Point vector = extractVectorFromAngle(arg);
                    if (!collides(player_cell, vector, nearby_cells, nearby_pheromes))
                        return new Move(vector, memory);
                }
                // if nothing worked, sit in place
                return new Move(new Point(0,0), memory);
            }
            
            return new Move(destination, memory);
        } else {
            /* Random walker strategy:
               Move away from friendly cells that are too close (specified by FRENDLY_AVOID_DIST)
               If no closeby friendly cells to avoid, act like default player (move in straight lines)
             */
            int prevAngle = memory & 0b01111111;

            Iterator<Cell> cell_it = nearby_cells.iterator();
            double sumX = 0;
            double sumY = 0;
            int count = 0;
            Point vector;

            // calculate avg position of nearby friendly cells
            while (cell_it.hasNext()) {
                Cell curr = cell_it.next();
                if (curr.player != player_cell.player)
                    continue;
                if (player_cell.distance(curr) > FRIENDLY_AVOID_DIST) // don't worry about far away friendlies
                    continue;
                Point currPos = curr.getPosition();
                sumX += currPos.x;
                sumY += currPos.y;
                count++;
            }

            if (count == 0) { // case: no friendly cells to move away from
                // if had a previous direction, keep going in that direction
                if (prevAngle > 0) {
                    vector = extractVectorFromAngle( (int)prevAngle);
                    if (!collides( player_cell, vector, nearby_cells, nearby_pheromes))
                        return new Move(vector, memory);
                }

                // if will collide or didn't have a previous direction, pick a random direction
                prevAngle = gen.nextInt(120);
                vector = extractVectorFromAngle(prevAngle);
            } else { // case: friendly cells too close, move in opposite direction
                double avgX = sumX / ((double) count);
                double avgY = sumY / ((double) count);

                double towardsAvgX = avgX - currentPosition.x;
                double towardsAvgY = avgY - currentPosition.y;

                double distanceFromAvg = Math.hypot(towardsAvgX, towardsAvgY);
            
                double awayX = (-(towardsAvgX)/distanceFromAvg) * Cell.move_dist;
                double awayY = (-(towardsAvgY)/distanceFromAvg) * Cell.move_dist;

                // clear the previous vector bits
                prevAngle = 0;
                vector = new Point(awayX, awayY);
            }
            
            if (collides(player_cell, vector, nearby_cells, nearby_pheromes)) {
                int arg = gen.nextInt(120);
                // try 20 times to avoid collision
                for (int i = 0; i < 20; i++) {
                    vector = extractVectorFromAngle(arg);
                    if (!collides(player_cell, vector, nearby_cells, nearby_pheromes)) {
                        memory = (byte) (0b10000000 | arg);
                        return new Move(vector, memory);
                    }

                }
                // if still keeps colliding, stay in place
                return new Move(new Point(0,0), (byte) 0b10000000);
            }

            // didn't collide, recreate memory with new vector and move
            memory = (byte) (0b10000000 | prevAngle);
            return new Move(vector, memory);
        }
    }

    // check if moving player_cell by vector collides with any nearby cell or hostile pherome
    private boolean collides(Cell player_cell, Point vector, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
	Iterator<Cell> cell_it = nearby_cells.iterator();
	Point destination = player_cell.getPosition().move(vector);
	while (cell_it.hasNext()) {
	    Cell other = cell_it.next();
	    if ( destination.distance(other.getPosition()) < 0.5*player_cell.getDiameter() + 0.5*other.getDiameter() + 0.00011)
		return true;
	}
	Iterator<Pherome> pherome_it = nearby_pheromes.iterator();
	while (pherome_it.hasNext()) {
	    Pherome other = pherome_it.next();
            if (other.player != player_cell.player && destination.distance(other.getPosition()) < 0.5*player_cell.getDiameter() + 0.0001)
                return true;
	}
	return false;
    }
    
    // convert an angle (in 3-deg increments) to a vector with magnitude Cell.move_dist (max allowed movement distance)
    private Point extractVectorFromAngle(int arg) {
	double theta = Math.toRadians( 3* (double)arg );
	double dx = Cell.move_dist * Math.cos(theta);
	double dy = Cell.move_dist * Math.sin(theta);
	return new Point(dx, dy);
    }

    private int getDirection(byte mem) {
            return (mem & 3);
    }

    private int getDuration(byte mem) {
            return ((mem >> 2) & 3);
    }
    private byte writeDirection(int direction, byte memory) {
            int actualDirection = direction % NUMDIRECTIONS;
            byte mem = (byte)((memory & 0b11111100) | actualDirection);
            return mem;
    }
    private byte writeDuration(int duration, byte memory, int maxDuration) {
            int actualDuration = duration % maxDuration;
            byte mem = (byte)((memory & 0b11110011) | (actualDuration << 2));
            return mem;
    }

    private int getMaxDuration(int t, int numdirs) {
            return Math.min((t / numdirs), DURATION_CAP);
    }

    private Point getNewDest(int direction) {

            if (direction == 0) {
                    return new Point(0*Cell.move_dist,-1*Cell.move_dist);

            } else if (direction == 1) {
                    return new Point(1*Cell.move_dist,0*Cell.move_dist);

            } else if (direction == 2) {
                    return new Point(0*Cell.move_dist,1*Cell.move_dist);

            } else if (direction == 3) {
                    return new Point(-1*Cell.move_dist,0*Cell.move_dist);

            } else {
                    return new Point(0,0);
            }

    }

    private int getStrategy(byte memory) {
        int strategy = (memory >> 7) & 1;
        return strategy;
    }
}
