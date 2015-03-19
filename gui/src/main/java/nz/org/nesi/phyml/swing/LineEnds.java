

    /*
     * Adjusts line endings.
     * Copyright (C) 2001-2010 Stephen Ostermiller
     * http://ostermiller.org/contact.pl?regarding=Java+Utilities
     *
     * This program is free software; you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation; either version 2 of the License, or
     * (at your option) any later version.
     *
     * This program is distributed in the hope that it will be useful,
     * but WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     * GNU General Public License for more details.
     *
     * See LICENSE.txt for details.
     */

    package nz.org.nesi.phyml.swing;

    import java.io.*;
    import java.text.MessageFormat;
    import java.util.ResourceBundle;
    import java.util.Locale;

    /**
     * Stream editor to alter the line separators on text to match
     * that of a given platform.
     * More information about this class is available from <a target="_top" href=
     * "http://ostermiller.org/utils/LineEnds.html">ostermiller.org</a>.
     *
     * @author Stephen Ostermiller http://ostermiller.org/contact.pl?regarding=Java+Utilities
     * @since ostermillerutils 1.00.00
     */
    public class LineEnds {

    	/**
    	 * Version number of this program
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public static final String version = "1.2";


    	/**
    	 * The system line ending as determined
    	 * by System.getProperty("line.separator")
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public final static int STYLE_SYSTEM = 0;
    	/**
    	 * The Windows and DOS line ending ("\r\n")
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public final static int STYLE_WINDOWS = 1;
    	/**
    	 * The Windows and DOS line ending ("\r\n")
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public final static int STYLE_DOS = 1;
    	/**
    	 * The Windows and DOS line ending ("\r\n")
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public final static int STYLE_RN = 1;
    	/**
    	 * The UNIX and Java line ending ("\n")
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public final static int STYLE_UNIX = 2;
    	/**
    	 * The UNIX and Java line ending ("\n")
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public final static int STYLE_N = 2;
    	/**
    	 * The UNIX and Java line ending ("\n")
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public final static int STYLE_JAVA = 2;
    	/**
    	 * The Macintosh line ending ("\r")
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public final static int STYLE_MAC = 3;
    	/**
    	 * The Macintosh line ending ("\r")
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public final static int STYLE_R = 3;

    	/**
    	 * Buffer size when reading from input stream.
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	private final static int BUFFER_SIZE = 1024;
    	private final static int STATE_INIT = 0;
    	private final static int STATE_R = 1;

    	private final static int MASK_N = 0x01;
    	private final static int MASK_R = 0x02;
    	private final static int MASK_RN = 0x04;

    	/**
    	 * Change the line endings of the text on the input stream and write
    	 * it to the output stream.
    	 *
    	 * The current system's line separator is used.
    	 *
    	 * @param in stream that contains the text which needs line number conversion.
    	 * @param out stream where converted text is written.
    	 * @return true if the output was modified from the input, false if it is exactly the same
    	 * @throws BinaryDataException if non-text data is encountered.
    	 * @throws IOException if an input or output error occurs.
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public static boolean convert(InputStream in, OutputStream out) throws IOException {
    		return convert(in, out, STYLE_SYSTEM, true);
    	}

    	/**
    	 * Change the line endings of the text on the input stream and write
    	 * it to the output stream.
    	 *
    	 * @param in stream that contains the text which needs line number conversion.
    	 * @param out stream where converted text is written.
    	 * @param style line separator style.
    	 * @return true if the output was modified from the input, false if it is exactly the same
    	 * @throws BinaryDataException if non-text data is encountered.
    	 * @throws IOException if an input or output error occurs.
    	 * @throws IllegalArgumentException if an unknown style is requested.
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public static boolean convert(InputStream in, OutputStream out, int style) throws IOException {
    		return convert(in, out, style, true);
    	}

    	/**
    	 * Change the line endings of the text on the input stream and write
    	 * it to the output stream.
    	 *
    	 * The current system's line separator is used.
    	 *
    	 * @param in stream that contains the text which needs line number conversion.
    	 * @param out stream where converted text is written.
    	 * @param binaryException throw an exception and abort the operation if binary data is encountered and binaryExcepion is false.
    	 * @return true if the output was modified from the input, false if it is exactly the same
    	 * @throws BinaryDataException if non-text data is encountered.
    	 * @throws IOException if an input or output error occurs.
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public static boolean convert(InputStream in, OutputStream out, boolean binaryException) throws IOException {
    		return convert(in, out, STYLE_SYSTEM, binaryException);
    	}

    	/**
    	 * Change the line endings of the text on the input stream and write
    	 * it to the output stream.
    	 *
    	 * @param in stream that contains the text which needs line number conversion.
    	 * @param out stream where converted text is written.
    	 * @param style line separator style.
    	 * @param binaryException throw an exception and abort the operation if binary data is encountered and binaryExcepion is false.
    	 * @return true if the output was modified from the input, false if it is exactly the same
    	 * @throws BinaryDataException if non-text data is encountered.
    	 * @throws IOException if an input or output error occurs.
    	 * @throws IllegalArgumentException if an unknown style is requested.
    	 *
    	 * @since ostermillerutils 1.00.00
    	 */
    	public static boolean convert(InputStream in, OutputStream out, int style, boolean binaryException) throws IOException {
    		byte[] lineEnding;
    		switch (style) {
    			case STYLE_SYSTEM: {
    				 lineEnding = System.getProperty("line.separator").getBytes();
    			} break;
    			case STYLE_RN: {
    				 lineEnding = new byte[]{(byte)'\r',(byte)'\n'};
    			} break;
    			case STYLE_R: {
    				 lineEnding = new byte[]{(byte)'\r'};
    			} break;
    			case STYLE_N: {
    				 lineEnding = new byte[]{(byte)'\n'};
    			} break;
    			default: {
    				throw new IllegalArgumentException("Unknown line break style: " + style);
    			}
    		}
    		byte[] buffer = new byte[BUFFER_SIZE];
    		int read;
    		int state = STATE_INIT;
    		int seen = 0x00;
    		while((read = in.read(buffer)) != -1){
    			for (int i=0; i<read; i++){
    				byte b = buffer[i];
    				if (state==STATE_R){
    					if(b!='\n'){
    						out.write(lineEnding);
    						seen |= MASK_R;
    					}
    				}
    				if (b=='\r'){
    					state = STATE_R;
    				} else {
    					if (b=='\n'){
    						if (state==STATE_R){
    							seen |= MASK_RN;
    						} else {
    							seen |= MASK_N;
    						}
    						out.write(lineEnding);
    					} else if(binaryException && b!='\t' && b!='\f' && (b & 0xff)<32){
    						throw new RuntimeException("Binary data encountered, line break replacement aborted.");
    					} else {
    						out.write(b);
    					}
    					state = STATE_INIT;
    				}
    			}
    		}
    		if (state==STATE_R){
    			out.write(lineEnding);
    			seen |= MASK_R;
    		}
    		if (lineEnding.length==2 && lineEnding[0]=='\r' && lineEnding[1]=='\n'){
    			return ((seen & ~MASK_RN)!=0);
    		} else if (lineEnding.length==1 && lineEnding[0]=='\r'){
    			return ((seen & ~MASK_R)!=0);
    		} else if (lineEnding.length==1 && lineEnding[0]=='\n'){
    			return ((seen & ~MASK_N)!=0);
    		} else {
    			return true;
    		}
    	}





    }


