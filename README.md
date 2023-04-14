# file-compressor
- My program implements file compression and decompression for UTF-8 format files using Huffman coding. The program reads in the input file, builds a Huffman tree based on the frequency of characters, and generates a variable-length encoding for each character. The encoded data is then written to an output file, which can be later decoded using the same Huffman tree.

- During the decompression phase, the program reads in the compressed file and the Huffman tree used to encode it, and then reads in the encoded data bit by bit, decoding it using the Huffman tree until the original data is recovered. The resulting decompressed data is then written to an output file.

- The use of Huffman coding allows for efficient compression of the data while minimizing data loss.
