package slather.g7;

// TODO This is very inefficient. Need to use bit operations to improve performance
public class Memory {
	public int opposite;
	public int moveDirectionCountdown;
	public int offsetCountDown;
	public int moveDirectionCountSetter;

	public Memory(byte memory) {
		String memStr = memoryToString(memory);
		this.opposite = getMemoryBlock(memStr, 7, 8);
		this.moveDirectionCountdown = getMemoryBlock(memStr, 4, 7);
		this.offsetCountDown = getMemoryBlock(memStr, 2, 4);	
		this.moveDirectionCountSetter = getMemoryBlock(memStr, 0, 2); // Initializes to 0 to start off
	}

	Memory(int offset, int dirCnt, int opposite) {
		this.offsetCountDown = offset;
		this.moveDirectionCountdown = dirCnt;
		this.opposite = opposite;
		this.moveDirectionCountSetter = 0;
	}

	public byte getByte() {
		String memStr = blockToString(moveDirectionCountSetter, 2) + 
						blockToString(offsetCountDown, 2) + 
						blockToString(moveDirectionCountdown, 3) +
						blockToString(opposite, 1);
		
		if (memStr.length() != 8) {
			System.out.println("The memory is converted to " + memStr.length() + " bits: " + memStr);
			return Byte.parseByte("0", 2);
		}
		
		int integer = Integer.parseInt(memStr, 2);
		if (integer > 127) {
			integer = 256 - integer;
		}
		byte b = (byte) integer;
		return b;
	}

	public static String memoryToString(byte memory) {
		String s = String.format("%8s", Integer.toBinaryString(memory & 0xFF)).replace(' ', '0');
		return s;
	}

	public static char getMemoryAt(String memoryStr, int index) {
		if (index > memoryStr.length())
			return ' ';
		return memoryStr.charAt(index - 1);
	}

	/*
	 * Retrieving a block of memory from index start inclusively to end
	 * exclusively
	 */
	public static int getMemoryBlock(String memoryStr, int start, int end) {
		String sub = memoryStr.substring(start, end);
		int val = Integer.parseInt(sub, 2);
		return val;
	}

	public static String blockToString(int block, int length) {
		String target = Integer.toBinaryString(block);
		if (target.length() > length) {
			System.out.println("The information " + block + " is too big to fit in memory of length " + length);
			return target.substring(0, length);
		} else if (target.length() < length) {
			StringBuilder sb = new StringBuilder(target);
			for (int i = target.length(); i < length; i++) {
				sb.insert(0, '0');
			}
			return sb.toString();
		} else {
			return target;
		}
	}
}