#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <random>
//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
namespace {
//---------------------------------------------------------------------------
// A sparse distance matrix, must be symmetric
class Matrix {
   public:
   /// An entry
   struct Entry {
      /// The column
      unsigned column;
      /// The weight
      double weight;
   };
   /// An entry range
   struct EntryRange {
      /// The entry range
      const Entry *from, *to;

      /// Empty range?
      bool empty() const { return from == to; }
      /// First element
      const Entry* begin() { return from; }
      /// Behind the last element
      const Entry* end() { return to; }
   };

   public:
   /// The entries
   vector<Entry> entries;
   /// The entries offsets
   vector<unsigned> entryOffsets;
   /// The width of the matrix
   unsigned width = 0;

   public:
   /// Get the number of rows
   unsigned getRowCount() const { return width; }
   /// Get all entries in a row
   EntryRange getRow(unsigned i) const { return EntryRange{entries.data() + entryOffsets[i - 1], entries.data() + entryOffsets[i]}; }

   /// Read a file
   static Matrix readFile(string fileName);
};
//---------------------------------------------------------------------------
Matrix Matrix::readFile(string fileName)
// Read a sparse matrix file
{
   ifstream in(fileName);
   if (!in.is_open()) {
      cerr << "unable to read " << fileName << endl;
      exit(1);
   }

   Matrix result;

   // Check the header line
   string line;
   while (getline(in, line)) {
      if (line.empty() || (line.front() == '%')) continue;

      stringstream s(line);
      unsigned width, height, entries;
      s >> width >> height >> entries;
      if ((width != height) || (!width) || (!entries)) {
         cerr << "matrix must be symmetric and non-empty" << endl;
         exit(1);
      }
      result.width = width;
      result.entries.reserve(entries);
      result.entryOffsets.reserve(width + 1);
      break;
   }

   // Read the elements
   unsigned currentRow = 0, currentColumn = 0;
   while (getline(in, line)) {
      if (line.empty() || (line.front() == '%')) continue;
      stringstream s(line);
      unsigned row, column;
      double weight;
      if (!(s >> row >> column >> weight)) {
         cerr << "malformed matrix line " << line << endl;
         exit(1);
      }
      if ((row > result.width) || (column > result.width)) {
         cerr << "malformed matrix format, cell offset out of bounds " << line << endl;
         exit(1);
      }
      if (row < currentRow) {
         cerr << "malformed matrix format, row number decreased " << line << endl;
         exit(1);
      }
      if (row == currentRow) {
         if (column <= currentColumn) {
            cerr << "malformed matrix format, column number decreased " << line << endl;
            exit(1);
         }
         currentColumn = column;
      } else {
         result.entryOffsets.insert(result.entryOffsets.end(), row - currentRow, result.entries.size());
         currentRow = row;
         currentColumn = column;
      }
      result.entries.push_back({column, weight});
   }
   result.entryOffsets.insert(result.entryOffsets.end(), result.width + 1 - currentRow, result.entries.size());
   return result;
}
//---------------------------------------------------------------------------
/// A priority queue of entries. We need both a min heap (best candidate) and a max heap (worst candidate), thus we implement it only once and use a template paramter
template <bool isMin>
class PriorityQueue {
   public:
   /// An entry
   struct Entry {
      /// The id
      unsigned index;
      /// The weight
      double weight;
   };

   public:
   /// The entries
   vector<Entry> entries;
   /// The entry lookup
   unordered_map<unsigned, unsigned> entryLookup;

   /// Comparison logic
   static bool isLess(double a, double b) {
      if (isMin) {
         return a < b;
      } else {
         return a > b;
      }
   }

   /// Move an element up in the heap until it is at the correct position
   void heapifyUp(unsigned slot);
   /// Move an element down the heap until it is at the correct position
   void heapifyDown(unsigned slot);

   public:
   /// Is the heap empty?
   bool empty() const { return entries.empty(); }
   /// The number of entries in the heap
   unsigned size() const { return entries.size(); }

   /// Get the top element
   const Entry& front() const { return entries.front(); }
   /// Remove the top element from the queue
   Entry pop_front();
   /// Add an element
   void insert(Entry entry);

   using iterator = const Entry*;
   /// Find an entry
   iterator find(unsigned index) const;
   /// Marker for not found
   iterator end() const { return nullptr; }
   /// Update an entry
   void update(iterator pos, double newWeight);
   /// Remove an entry
   void erase(iterator pos);
};
//---------------------------------------------------------------------------
template <bool isMin>
void PriorityQueue<isMin>::heapifyUp(unsigned slot)
// Move an element up in the heap until it is at the correct position
{
   if ((!slot) || (slot >= entries.size())) return;
   auto& currentPos = entryLookup[entries[slot].index];

   // Bubble up until the heap condition is restored
   while (slot > 0) {
      unsigned parentSlot = slot / 2;
      if (isLess(entries[slot].weight, entries[parentSlot].weight)) {
         entryLookup[entries[parentSlot].index] = slot;
         currentPos = parentSlot;
         swap(entries[slot], entries[parentSlot]);
         slot = parentSlot;
      } else {
         break;
      }
   }
}
//---------------------------------------------------------------------------
template <bool isMin>
void PriorityQueue<isMin>::heapifyDown(unsigned slot)
// Move an element down the heap until it is at the correct position
{
   if (slot >= entries.size()) return;
   auto& currentPos = entryLookup[entries[slot].index];

   // Bubble down until the heap condition is restored
   while (true) {
      unsigned leftChild = 2 * slot, rightChild = leftChild + 1;
      unsigned selectedChild;
      if (rightChild < entries.size()) {
         selectedChild = isLess(entries[leftChild].weight, entries[rightChild].weight) ? leftChild : rightChild;
      } else if (leftChild < entries.size()) {
         selectedChild = leftChild;
      } else {
         break;
      }
      if (isLess(entries[selectedChild].weight, entries[slot].weight)) {
         entryLookup[entries[selectedChild].index] = slot;
         currentPos = selectedChild;
         swap(entries[slot], entries[selectedChild]);
         slot = selectedChild;
      } else {
         break;
      }
   }
}
//---------------------------------------------------------------------------
template <bool isMin>
typename PriorityQueue<isMin>::Entry PriorityQueue<isMin>::pop_front()
// Remove the top element from the queue
{
   auto result = entries.front();
   swap(entries.front(), entries.back());
   entryLookup[entries.front().index] = 0;
   entryLookup.erase(result.index);
   entries.pop_back();
   heapifyDown(0);

   return result;
}
//---------------------------------------------------------------------------
template <bool isMin>
void PriorityQueue<isMin>::insert(Entry entry)
// Add an element
{
   unsigned slot = entries.size();
   entries.push_back(entry);
   entryLookup[entry.index] = slot;
   heapifyUp(slot);
}
//---------------------------------------------------------------------------
template <bool isMin>
typename PriorityQueue<isMin>::iterator PriorityQueue<isMin>::find(unsigned index) const
// Find an entry
{
   auto iter = entryLookup.find(index);
   if (iter == entryLookup.end()) return nullptr;
   return entries.data() + iter->second;
}
//---------------------------------------------------------------------------
template <bool isMin>
void PriorityQueue<isMin>::update(iterator pos, double newWeight)
// Update an entry
{
   unsigned slot = pos - entries.data();
   if (isLess(newWeight, pos->weight)) {
      entries[slot].weight = newWeight;
      heapifyUp(slot);
   } else if (isLess(pos->weight, newWeight)) {
      entries[slot].weight = newWeight;
      heapifyDown(slot);
   }
}
//---------------------------------------------------------------------------
template <bool isMin>
void PriorityQueue<isMin>::erase(iterator pos)
// Remove an entry
{
   unsigned index = pos->index;
   unsigned slot = pos - entries.data();
   double oldWeight = pos->weight, newWeight=entries.back().weight;
   swap(entries[slot], entries.back());
   entryLookup[entries[slot].index] = slot;
   entryLookup.erase(index);
   entries.pop_back();
   if (isLess(oldWeight,newWeight)) {
      heapifyDown(slot);
   } else {
      heapifyUp(slot);
   }
}
//---------------------------------------------------------------------------
using MinHeap = PriorityQueue<true>;
using MaxHeap = PriorityQueue<false>;
//---------------------------------------------------------------------------
static vector<Matrix::Entry> getTopKNeighbors(const Matrix& m, unsigned start, unsigned k)
// Find the top k neighbors for a node
{
   vector<Matrix::Entry> result;
   result.reserve(k);

   MinHeap minHeap;
   MaxHeap maxHeap;
   minHeap.insert({start, 0});
   maxHeap.insert({start, 0});
   while (!minHeap.empty()) {
      // Remove the next element and add it to result
      auto current = minHeap.pop_front();
      if (current.index != start) {
         result.push_back({current.index, current.weight});
         if (result.size() >= k) break;
      }

      // Examine all outgoing edges
      for (auto& e : m.getRow(current.index)) {
         // Already in heap?
         double newDistance = current.weight + e.weight;
         auto iter = maxHeap.find(e.column);
         if (iter != maxHeap.end()) {
            // If an entry is in the max heap but not in the min heap it is already in the result and we can ignore it
            auto iter2 = minHeap.find(e.column);
            if (iter2 != minHeap.end()) {
               // Update if shorter
               if (newDistance < iter2->weight) {
                  minHeap.update(iter2, newDistance);
                  maxHeap.update(iter, newDistance);
               }
            }
         } else if (maxHeap.size() <= k) {
            // As long as we have not seen k candidates add it unconditionally
            minHeap.insert({e.column, newDistance});
            maxHeap.insert({e.column, newDistance});
         } else if (newDistance < maxHeap.front().weight) {
            // We got a better candidate, remove the old one
            auto worst = maxHeap.pop_front();
            minHeap.erase(minHeap.find(worst.index));
            minHeap.insert({e.column, newDistance});
            maxHeap.insert({e.column, newDistance});
         }
      }
   }

   return result;
}
//---------------------------------------------------------------------------
}
//---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
   if (argc != 6) {
      cout << "usage: " << argv[0] << " matrixfile, output_prefix, k, n_chunks, percent_subsample" << endl;
      return 0;
   }

   // Read the matrix file
   Matrix matrix = Matrix::readFile(argv[1]);
   //get output_prefix
   string output_prefix = argv[2];
   // The number of neighbors we are interested in
   unsigned k = stoi(argv[3]); //convert input char to integer

   unsigned n_chunks = stoi(argv[4]);
   unsigned limit = matrix.getRowCount();
   unsigned len_ch;
   if (n_chunks <= 1){
	   n_chunks = 1;
	   len_ch = limit;
   }
   else{
   	len_ch = limit / (n_chunks - 1);
   }
   //get percentage to which should be subsampled
   unsigned sub = stoi(argv[5]);

   //ininitialize random number generator
   random_device rd; //used to obtain seed for random number engine
   mt19937 gen(rd()); //standard merseen_twister_engine seeded with rd()
   uniform_int_distribution<> dis(1, 100); //uniform int distribution between 0 and 100
   int rand_res;
   //sanity check
   //double sum = 0;

   //variable declaration
   ofstream distances;
   ofstream indices;
   string dist;
   string indi;
   unsigned lower;
   unsigned upper;

   for (unsigned n_ch = 0; n_ch < n_chunks; ++n_ch){
   // Find the top k elements for all nodes. Computes the sum of all the weights, just to have some result to show
   // write all neighbors and weights to two files
    dist = output_prefix + "_distances_" + to_string(n_ch) + ".txt";
    indi = output_prefix + "_indices_" + to_string(n_ch) + ".txt";

    distances.open(dist, ios::out | ios::binary);
    indices.open(indi, ios::out | ios::binary);
    lower = n_ch * len_ch + 1;
    upper = (n_ch+1)*len_ch;
    // don't run over upper limit
    if (upper > limit) {
     upper = limit;
    }
    for (unsigned row = lower; row <= upper; row++) {
      //cout << row << endl;
      // Ignore empty rows
      if (matrix.getRow(row).empty()) {
	  //distances << row << endl; //add index to distances file to keep order
	  //indices << row << endl; //add index to indices file to keep order
	  continue;
      }
      // use subsampling
      rand_res = dis(gen); //generate random number
      if (rand_res> sub){ //skip for 100-sub percent of the data
          continue;
      }

      // Find the top k neighbors
      auto neighbors = getTopKNeighbors(matrix, row, k);
      distances << row; //add index of the root index to file
      indices << row;   //add index of the root index to file
      for (auto& n : neighbors) {
         distances << ',' << n.weight;
	 indices << ',' << n.column;
	 //sum += n.weight;
      }
      distances << endl; //add end of line after each knn search
      indices << endl;  //add end of line after each knn search
    }
    distances.close();
    indices.close();
   }
   //cout << "sum of weights for all top " << k << " neighbors " << sum << endl;
}
//---------------------------------------------------------------------------
