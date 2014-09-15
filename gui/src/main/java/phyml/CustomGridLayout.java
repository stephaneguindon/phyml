package phyml;

import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.LayoutManager;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;

/**
 * This LayoutManager specifies the width and height of components by fractions
 * specified by the user. For example the commands<br>
 * <br>
 * 
 * CustomGridLayout layout = new CustomGridLayout();<br>
 * setLayout(layout);<br>
 * layout.setDimensions(0.5,0.5);<br>
 * add(comp)<br>
 * <br>
 * 
 * specifies the width and height of the the component as half the width and
 * half the height of the parent container. The next component would be added at
 * the upper right corner of former if the specified width and height allows
 * that. A component which is to big to be added is not added. Therefore the
 * total width and heights specified in the setDimensions() method need to be
 * smaller or equal 1. <br>
 * 
 * @author Christoph Knapp
 * @version 15-June-2012
 */
public class CustomGridLayout implements LayoutManager {
	private ArrayList<Rectangle> orientations;
	private ArrayList<Point> possPoints;
	private ArrayList<Area> areas;
	private final double HORIZONTAL = 0.2;
	private final double VERTICAL = 0.2;

	/**
	 * Constructor method for CustomGridLayout.
	 */
	public CustomGridLayout() {
		orientations = new ArrayList<Rectangle>();
		possPoints = new ArrayList<Point>();
		areas = new ArrayList<Area>();
	}

	@Override
	public void layoutContainer(Container parent) {
		// TODO Auto-generated method stub
		setup(parent);
		int ncomponents = parent.getComponentCount();
		for (int i = 0; i < ncomponents; i++) {
			Component comp = parent.getComponent(i);
			comp.setBounds(orientations.get(i));
		}
	}

	/**
	 * Sets up the layout components. There start locations and sizes.
	 * 
	 * @param parent
	 *            Container : The parent container of the component added.
	 */
	private void setup(Container parent) {
		possPoints = new ArrayList<Point>();
		possPoints.add(new Point(0, 0));
		orientations = new ArrayList<Rectangle>();
		int ncomponents = parent.getComponentCount();
		if (areas.size() < ncomponents) {
			for (int i = areas.size(); i < ncomponents; i++) {
				if (areas.size() > 0) {
					areas.add(new Area(areas.get(areas.size() - 1)
							.getHorizontal(), areas.get(areas.size() - 1)
							.getVertical()));
				} else {
					areas.add(new Area(HORIZONTAL, VERTICAL));
				}
			}
		}
		for (int i = 0; i < ncomponents; ++i) {
			Rectangle r = new Rectangle(0, 0, (int) (parent.getWidth() * areas
					.get(i).getHorizontal()), (int) (parent.getHeight() * areas
					.get(i).getVertical()));
			boolean added = false;
			Point p = null;
			for (int j = 0; j < possPoints.size(); j++) {
				p = possPoints.get(j);
				r.setLocation(p);
				boolean hasIntersection = false;
				for (int k = 0; k < orientations.size(); k++) {
					Rectangle add = orientations.get(k);
					if (r.intersects(add)) {
						hasIntersection = true;
					}
				}
				if (hasIntersection) {
					continue;
				}
				if (p.x + r.width <= parent.getWidth()
						&& p.y + r.height <= parent.getHeight()) {
					orientations.add(i, r);
					added = true;
					break;
				}
			}
			if (!added) {
				orientations.add(i, new Rectangle(0, 0, 0, 0));
			} else {
				if (p.x + r.width < parent.getWidth()) {
					possPoints.add(new Point(p.x + r.width, p.y));
				}
				if (p.y + r.height < parent.getHeight()) {
					possPoints.add(new Point(p.x, p.y + r.height));
				}
				order();
				possPoints.remove(p);
			}
		}
	}

	/**
	 * Simple divide and conquer sorting algorythm for Point objects. The points
	 * need to be ordered or things become added to the wrong positions.
	 */
	private void order() {
		ArrayList<Point> temp = new ArrayList<Point>();
		for (int i = 0; i < possPoints.size(); i++) {
			if (i == 0) {
				temp.add(possPoints.get(i));
			} else {
				int count = 0;
				boolean added = false;
				for (Point p : temp) {
					if (compare(p, possPoints.get(i)) == 1) {
						temp.add(count, possPoints.get(i));
						added = true;
						break;
					} else if (compare(p, possPoints.get(i)) == 0) {
						break;
					}
					count++;
				}
				if (!added) {
					temp.add(possPoints.get(i));
				}
			}
		}
		possPoints = temp;
	}

	/**
	 * Part of the order() method to compare two Point objects.
	 * 
	 * @param point1
	 *            Point : Point 1 to be compared.
	 * @param point2
	 *            Point : Point 2 to be compared.
	 * @return int : -1 if point1 is smaller than point2, 0 if they are the same
	 *         and 1 if point1 is bigger than point2.
	 */
	private int compare(Point point1, Point point2) {
		if (point1.y < point2.y) {
			return -1;
		} else if (point1.y > point2.y) {
			return 1;
		} else {
			if (point1.x < point2.x) {
				return -1;
			} else if (point1.x > point2.x) {
				return 1;
			}
		}
		return 0;
	}

	@Override
	public void addLayoutComponent(String name, Component comp) {
		// TODO Auto-generated method stub
		System.out.println(name + comp.toString() + "addLayoutComponent");
	}

	@Override
	public Dimension minimumLayoutSize(Container parent) {
		// TODO Auto-generated method stub
		System.out.println(parent.toString() + "minimumLayoutSize");
		return null;
	}

	@Override
	public Dimension preferredLayoutSize(Container parent) {
		// TODO Auto-generated method stub
		System.out.println(parent.toString() + "preferredLayoutSize");
		return null;
	}

	@Override
	public void removeLayoutComponent(Component comp) {
		// TODO Auto-generated method stub
		System.out.println(comp.toString() + "removeLayoutComponent");
	}

	/**
	 * This method sets the dimension of a component inside its parent
	 * container. If the dimension is to big to add another component, the
	 * component is not added.<br>
	 * 
	 * @param horizontal
	 *            double value <=1
	 * @param vertical
	 *            double value <=1
	 */
	public void setDimensions(double horizontal, double vertical) {
		if (horizontal <= 1 && vertical <= 1) {
			areas.add(new Area(horizontal, vertical));
		}
	}

	@Override
	public String toString() {
		String s = "";
		for (Area a : this.areas) {
			s = s + "\n" + a.toString();
		}
		for (Rectangle r : this.orientations) {
			s = s + "\n" + r.toString();
		}
		for (Point p : this.possPoints) {
			s = s + "\n" + p.toString();
		}
		return s;
	}

	/**
	 * Class for specifying the horizontal and vertical proportion of the
	 * component in the gui.
	 * 
	 * @author Christoph Knapp
	 * 
	 */
	private class Area {
		private double horizontal;
		private double vertical;

		/**
		 * Constructor method
		 * 
		 * @param horizontal
		 *            double : horizontal proportion of the added component in
		 *            the parent container. The value of it should be >0 and
		 *            <=1.
		 * @param vertical
		 *            double : vertical proportion of the added component in the
		 *            parent container. The value of it should be >0 and <=1.
		 */
		public Area(double horizontal, double vertical) {
			this.horizontal = horizontal;
			this.vertical = vertical;
		}

		/**
		 * Retrieves the horizontal proportion of a component.
		 * 
		 * @return double : horizontal
		 */
		public double getHorizontal() {
			return horizontal;
		}

		/**
		 * Retrieves the vertical proportion of a component.
		 * 
		 * @return double : vertical
		 */
		public double getVertical() {
			return vertical;
		}

		@Override
		public String toString() {
			return "Area [horizontal=" + horizontal + ", vertical=" + vertical
					+ "]";
		}
	}
}
