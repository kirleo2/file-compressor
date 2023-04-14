#ifndef __PROGTEST__
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <climits>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <set>
#include <queue>
#include <memory>
#include <functional>
#include <stdexcept>
using namespace std;



// Note: classes and methods are commented outsite, functions are commented inside.

/*
 First half of the code is responsible for decompression.
 */



/*
 Class for reading file bits and bytes from input file
 */
class CBitReader {
    
public:
    /*
     reading input file bytes, those we will keep in vector of unsigned 8-bit integers
     (the same thing as an unsigned char)
     */
    CBitReader ( const char * filename ) {
   
        ifstream file;
        file.open(filename, ios::binary);
        file_is_open = file.is_open();
        char byte;
        while (file.get(byte)) {
            uint8_t b = (int) (unsigned char) byte;
            file_bytes.push_back(b);
        }
        file.close();
    }
    bool isFileEmpty() {
        return !file_bytes.size();
    }
    /* examination, if we can proceed with getting of next bytes and bits in our input file
     (if we didnt reach the end of the file)
     */
    bool isBytesLeft() const {
        return (current_byte_pos < (int) file_bytes.size());
    }
    
    bool isFileOpen () const {
        return file_is_open;
    }
    
    /*
     basis of the program, reading bits using the mask 0b10000000
     we move the "1" towards the right and get current bit in current byte using & operator
     */
    bool getNextBit ( ) {
       
        bool bit = ((file_bytes[current_byte_pos]) & ( mask >> (current_bit_pos++)));
        if (current_bit_pos == 8) {
            current_bit_pos = 0;
            current_byte_pos++;
        }
        return bit;
    }
    
    /*
     Assemble the byte using bits and powers of 2
     */
    uint8_t getNextByte ( ) {
        
        uint8_t current_byte = 0;
        for (int i = 7; i >= 0 ;i--) {
            current_byte += getNextBit() * (1 << i);
        }
        return current_byte;
    }
    
    /*
     this method is based on getByte method and validation of bytes according to the unicode
     if it is not succeed, will be returned empty vector
     */
    vector<uint8_t> getNextUtf8Char ( ) {
       
        vector<uint8_t> utf8;
        uint8_t fb = getNextByte();
        int bit_pos = 0;
        // to get utf-8 char, we need to count the quantity of "1" in first byte
        while (fb & (mask >> bit_pos)) {
            bit_pos++;
        }
        
        if (bit_pos <= 4 && bit_pos != 1) {
            
            // check, if we are not out of the unicode borders
            
            if (fb == 0xf4) {
                uint8_t next_byte = getNextByte();
                if ((next_byte & 0x90) == 0x90 || (next_byte & 0xa0) == 0xa0) {
                    utf8.clear();
                    return utf8;
                }
                current_byte_pos--;
            }
            utf8.push_back(fb);
        for (int i = 1; i < bit_pos; i++) {
            uint8_t next_byte = getNextByte();
            if (isValidUtf8Byte(next_byte)) utf8.push_back(next_byte);
            else {
                utf8.clear();
                break;
            }
        }
        }
        return utf8;
        
    }
    /*
     overloaded operator, needed only to check the files are same
     */
    bool operator==(const CBitReader & other) {
     
        return file_bytes == other.file_bytes;
    }

private:
    int current_bit_pos = 0;
    int current_byte_pos = 0;
    bool file_is_open = false;
    const uint8_t mask = 0b10000000;
    vector <uint8_t> file_bytes;
private:
    /*
     validation of 2-4 bytes, they must satisfy the format 10xxxxxx
     */
    bool isValidUtf8Byte (uint8_t byte) {
       
        if ((byte & mask) && !(byte & (mask >> 1)) ) return true;
        else return false;
    }
};

/*
 Class, that is responsible for tree node
 it has pointers to left and right child
 and as data contains UTF-8 char (it will be set only if the node is a leaf)
 */
class TNode {
    
public:
    TNode () {
        left = nullptr;
        right = nullptr;
    }
    TNode * left, * right;
    vector<uint8_t> value;
};

/*
 Class, that represents tree root.
 It has pointers to left and right child
 */
class TTree {
    
public:
    TTree () {
        left = nullptr;
        right = nullptr;
    }
    /* free allocated memory*/
    ~TTree () {
       
        deleteNode(left);
        deleteNode(right);
    }
    TNode * left;
    TNode * right;
    
    /*
     Method, that is responsible for building of the tree
     if it is not succeed, returns false
     --------------------------------------
     br - object for reading bits and bytes
     --------------------------------------
     */
    bool buildTree(CBitReader & br) {
        if (!br.isBytesLeft()) return false;
        bool bit = br.getNextBit();
        // first bit in file bit must be "0" according to the representation of the tree
        if (bit) return false;
        
        left = new TNode();
        right = new TNode();
        bool left_node = buildNode(left, br);
        bool right_node = buildNode(right, br);
        return left_node & right_node;
        
    }
    
private:
    
    /*
     Recursive method to build all tree nodes.
     --------------------------------------
     node - pointer to the allocated parent
     br - object for reading bits and bytes
     --------------------------------------
     */
    bool buildNode(TNode * node, CBitReader & br) {
       
        
        if (!br.isBytesLeft()) return false;
        bool bit = br.getNextBit();
        // if bit is "1" then make a tree leaf with UTF-8 char
        if (bit) {
            
            vector<uint8_t> val = br.getNextUtf8Char();
            if (!val.size()) return false;
            node->value = val;
  
            return true;
        }
        // if bit is "0", then call yourself to make left and right child
        else {
            node->left = new TNode();
            node->right = new TNode();
            
            return (buildNode(node->left, br)) & (buildNode(node->right, br));
        }
    }
    void deleteNode(TNode * node) {
        if (node) {
            deleteNode(node->left);
            deleteNode(node->right);
            delete node;
        }
    }
};

bool getChunk(vector<uint8_t> & output_bytes, const TTree * root, CBitReader & br, int size) {
    /*
     Function, that is responsible for reading the chunk
     (decode all sequences using built tree)
     In a cycle we move through a tree, until we meet a leaf with data
     --------------------------------------
     output_bytes - vector of bytes of decompressed file
     root - pointer to the tree root
     br - object for reading bits and bytes
     size - size of a chunk
     --------------------------------------
     */
    for (int i = 0; i < size; i++) {
        // we must control, if the end of a file is not met yet
        if (!br.isBytesLeft()) return false;
        bool bit = br.getNextBit();
        
        TNode * ptr = !bit ? root->left : root->right;
        // move until the child exists, it could be either right or left
        while(ptr->left) {
            if (!br.isBytesLeft()) return false;
            bit = br.getNextBit();
            ptr = !bit ? ptr->left : ptr->right;
        }
        // put our UTF character into output bytes vector
        for (uint8_t byte : ptr->value) output_bytes.push_back(byte);
     
    
    }
    return true;
    
}

bool decompressFile ( const char * inFileName, const char * outFileName )
{
    /*
     The main task of the progtest. We need to decompress compressed text file.
     At first we need to read and build the binary tree, that encodes all UTF characters.
     TTree and TNode represents the tree.
     CBitReader represents bit stream of input file.
     Then we need do read the chunks and decode every character using our tree and assemble
     decompressed file at the end.
     */
    
    CBitReader br = CBitReader(inFileName);
    // if file is empty, then return false
    if (br.isFileEmpty()){
        return false;
    }
    // if file openning is not succeed, then return false
    if (!br.isFileOpen()){
        return false;}
    // class for Tree root
    TTree * prefix_tree = new TTree();
    bool tree_is_build = prefix_tree->buildTree(br);
    if (!tree_is_build) {
        return false;
        
    }
    
    // vector for output file bytes
    vector<uint8_t> output_bytes;
    // reading chunks with size of 4096 characters, unless we met 0
    while(br.getNextBit()) {
        if(!getChunk(output_bytes, prefix_tree, br, 4096)) return false;
    }
    // if we met 0, then read an integer - size of next chunk
    int chunk_size = 0;
    for (int i = 11; i >= 0; i--) {
        if (!br.isBytesLeft()) return false;
        bool bit = br.getNextBit();
        // assemble an integer from bits using powers of 2, 1 is 2^0, in cycle we increase the power
        chunk_size += bit * (1 << i);
    }
    // reading last chunk
    if(!getChunk(output_bytes, prefix_tree, br, chunk_size)) return false;

    
    // making the decompressed file
    ofstream output_file;
    output_file.open(outFileName,ios::binary);
    if (!output_file.is_open()) return false;
    for (uint8_t byte : output_bytes) {
        char b = (char) byte;
        if(!output_file.put(b)) return false;
    }
    output_file.close();
    // free a memory allocated for a tree root and all nodes (destructor will be called automatically)
    delete prefix_tree;
    
  return true;
}

/*
 
 
 
 The next part of code is responsible for compression
 
 
 
*/



/*
 Class, that represents the node of a encodying tree
 It has pointers to left and right child, parent and quantity of parents, and frequancy (for encodying algorihtm)
 */
class CNode {

public:
    /*
     Constructor to create a tree node
     */
    CNode(CNode *left, CNode * right) :
    left(left),
    right(right) {
        frequency = 0;
        parent_count = 1;
        parent = nullptr;
    }
    /*
     Constructor to create a tree leaf
     */
    CNode(vector<uint8_t> utf, int frequency):
    utf_value(utf),
    frequency(frequency){
        parent_count = 1;
        left = nullptr;
        right = nullptr;
        parent = nullptr;
        
    }
    vector<uint8_t> utf_value;
    int frequency;
    int parent_count;
    CNode * parent;
    CNode * left;
    CNode * right;
    
    
};

/*
 Class - the opposite of CBitReader.
 Its methods recieve bits and bytes to make a compressed file.
 */
class FileMaker {
   
public:
    FileMaker() {
        bit_position = 0;
        byte_position = 0;
        file_bytes.push_back(0x0);
    }
    /*
     This method writes bit on corrent position in current byte array using mask.
     We can imagine an infinite array of bytes 0x0 (0b00000000), that method fills it up
     with "1" on correct positions by moving "1" towards the right in mask.
     */
    void putBit(bool bit) {
        
        if (bit) {
        file_bytes[byte_position] |= (mask >> bit_position);
        }
        bit_position++;
        if (bit_position == 8) {
            bit_position = 0;
            byte_position++;
            file_bytes.push_back(0);
            //cout << "   " << endl;
        }
    }
    /*
     This method based on putBit. We extract bits using mask from received byte
     and put it in output array.
     */
    void putByte(uint8_t byte) {
  
        for (int i = 0; i < 8; i++) {
            putBit(byte & (mask >> i));
        }
    }
    /*
     Based on putByte.
     */
    void putUTF8(vector<uint8_t> utf) {
        
        for (uint8_t byte : utf) {
            putByte(byte);
        }
    }
    /*
     This method receives all encodying tree leafs and input file data.
     For each one of UTF character from input file bytes it goes through the tree (from leaf to the root)
     and records the path sequence of bits (for each node it checks from where it goes). After that the
     sequence has to be inverted (we need to write the path from root towards the leaf).
     --------------------------------------
     leafs - vector of tree leafs with data (UTF characters)
     file - vector of original (not encoded) input file charactaers
     --------------------------------------
     */
    void makeChunk(vector<CNode *> & leafs, vector<vector <uint8_t> > & file) {
     
        
        int bytes_left = (int ) file.size();
        
        // count of characters, that we have already passed through
        int counter = 0;
        
        // mask to build 12-bit integer (0b1000000000000)
        int b_mask = 0x800;
        
        for (vector <uint8_t> utf_char: file) {
            /*
             If we have passed through multiple of 4096 characters, we have to check how much remains.
             If more than 4096 - put "1" before the chunk. If less - put "0" and make a binary (12-bit) form
             of an amount of left characters.
             */
            if (counter % 4096 == 0) {
                counter=0;
                if (bytes_left >= 4096) {
                    putBit(true);
                }
                else {
                    putBit(false);
                    for(int i = 0; i < 12; i++) {
                        putBit(bytes_left & (b_mask >> i));
                    }
                }
            }
            CNode * node = nullptr;
            CNode * tmp = nullptr;
            // buffer in that we write the path from leaf to the root
            string buffer = "";
            
            // We need to find corresponding UTF leaf value
            for (CNode * leaf : leafs) {
                if (leaf->utf_value == utf_char) {
                    node = leaf;
                    break;
                }
            }
            // while the parent exists, we are moving towards the root
            while (node->parent) {
                // we need to remember from where we came
                tmp = node;
                node = node -> parent;
                buffer.append(node -> left == tmp ? "0" : "1");
            }
            // Passing through the binary path from end to the start to put the bits into the file
            for (int i = (int) buffer.length() - 1; i >= 0; i--) {
                if (buffer[i] == '1') putBit(true);
                else putBit(false);
            }
            
            counter++;
            bytes_left--;
            
            
        }
        
    }
    
    /*
     Method to write all encoded bytes into the result compressed file.
     */
    bool makeOutputFile(const char * filename) {
       ofstream out;
        out.open(filename);
        if(!out.is_open()) return false;
        for (uint8_t byte : file_bytes) {
            if(!out.put((char) byte)) return false;
        }
        out.close();
        return true;
    }

private:
    int bit_position;
    int byte_position;
    vector<uint8_t> file_bytes;
    uint8_t mask = 0b10000000;
    
};


bool compare(const CNode * c1, const CNode * c2) {
    /*
     Function for sorting array of nodes to build an encodying tree.
     */
    if (c1->frequency > c2->frequency) {
        return true;
    }
    else if (c1->frequency == c2->frequency) {
        if (c1->parent_count > c2->parent_count) return true;
        else return false;
    }
    else return false;
}

void serialize_tree(CNode * node, bool lr, vector<uint8_t> & queue, FileMaker & maker){
    /*
     Function of encodying the tree with pre-order method.
     --------------------------------------
     node - pointer to current node where we are located now
     queue - is UTF8 char, that we get from the leaf.
     lr - from where (left or right) the function was called.
     maker - object that is responsible for creating the output compressed file
     --------------------------------------
     */
    maker.putBit(lr);
    if( lr) {
        maker.putUTF8(queue);
        queue.clear();
    }
    if (node->left) {
        serialize_tree(node->left, false, queue, maker);
        serialize_tree(node->right, true, queue, maker);
    }
    else {
        queue = (node->utf_value);
    }
}

bool compressFile ( const char * inFileName, const char * outFileName )
{
    /*
     The main idea of compression is what we want to encode the most common characters
     with the shortest sequences.
     We need to count all occurancies of each UTF char, and then build an encodying tree.
     */
    
    CBitReader br (inFileName);
    if(!br.isFileOpen()) return false;
    // vector of UTF-8 chars, we need to count the occurances each one of them
    vector <vector <uint8_t> > utf8_bytes;
    while(br.isBytesLeft()) {
        vector<uint8_t> h = br.getNextUtf8Char();
        if (!h.size()) return false;
        utf8_bytes.push_back(h);
    }
    
    if (!utf8_bytes.size()) return false;
    /*
     We need make a copy, because we want to count the occurancies and at the same time
     we need to have a original file data
     */
    vector <vector <uint8_t> > utf8_bytes_copy;
    utf8_bytes_copy = utf8_bytes;
    // We will sort all chars and count which are next to each other
    sort(utf8_bytes.begin(), utf8_bytes.end());
   
    
    vector<CNode *> avaliable_nodes;
    int counter = 1;
    utf8_bytes.push_back(vector<uint8_t> ());
    for (int i = 1; i < (int) utf8_bytes.size(); i++) {
        if (utf8_bytes[i] == utf8_bytes[i-1]) counter++;
        else {
            // Every unique character (with its frequency) we put as a possible node in queue
            avaliable_nodes.push_back(new CNode(utf8_bytes[i - 1], counter));
            counter = 1;
        }
    }
    // We must have at least 2 unique characters
    if (avaliable_nodes.size() == 1) return false;
    
   
    //sort(avaliable_nodes.begin(), avaliable_nodes.end(), compare);
    
 
    vector<CNode *> leafs = avaliable_nodes;
/*
 We chose 2 nodes with a shortest frequency and a fewest parent's count and join them.
 Therefore we need to sort a vector in descending order in each iteration, and take
 two last elements.
 The encodying tree will be built from leafs towards the root.
 */
    while(avaliable_nodes.size() > 1) {
        sort(avaliable_nodes.begin(), avaliable_nodes.end(), compare);

        CNode * node1 = avaliable_nodes.back();
        avaliable_nodes.pop_back();

        CNode * node2 = avaliable_nodes.back();
        avaliable_nodes.pop_back();

        CNode * new_node = new CNode(node1, node2);
        new_node->frequency = node1->frequency + node2->frequency;
        new_node->parent_count = node1->parent_count + node2->parent_count;
        node1->parent = new_node;
        node2->parent = new_node;
        avaliable_nodes.push_back(new_node);

    }
    // Remaining node is the tree root
    CNode * root = avaliable_nodes.back();
    
    FileMaker maker = FileMaker();
    vector<uint8_t> queue;
    // write the tree into the file
    serialize_tree(root->left, false, queue, maker);
    serialize_tree(root->right, true, queue, maker);
    // last "1" defines the end of the tree and last character after itself.
    maker.putBit(true);
    maker.putUTF8(queue);
    // make chunks from input bytes using our encodying tree
    maker.makeChunk(leafs, utf8_bytes_copy);
    // put resulting bytes into compressed file
    if(!maker.makeOutputFile(outFileName)) return false;

    
   return true;
}
bool identicalFiles ( const char * fileName1, const char * fileName2 )
{
    CBitReader f1 = CBitReader(fileName1);
    CBitReader f2 = CBitReader(fileName2);
    return f1 == f2;
}

int main ( void )

{
    return 0;
}
