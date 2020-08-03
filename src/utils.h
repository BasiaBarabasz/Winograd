#include <vector>

#define DEBUG 0

using namespace std;

typedef vector<float> row;
typedef vector<row> matrix;
typedef vector<double> rowd;
typedef vector<rowd> matrixd;

template<typename T>
using Row = vector<T>;
template<typename T>
using Matrix = vector<Row<T>>;

/*
 * Class TreeNode
 */
template<typename T>
class TreeNode {
  public:
    /// order of magnitude, used in priority queue
    double key;
    /// summation result in subtree
    T add_res;

    /// pointer to left subtree
    TreeNode *left;
    /// pointer to right subtree
    TreeNode *right;

    /*
     * TreeNode constructor
     *
     * @param k order of magnitude
     * @param f element stored in a node
     */
    TreeNode(double k, T f) : key(k), add_res(f), left(NULL), right(NULL) { }
    /// getter for key
    double element() { return key; }
    /// getter for add_res
    T result() { return add_res; }
};

/*
 * Function prints given tree in 'inorder' order
 *
 * @param p pointer to root of a tree to print
 */
template<typename T>
void inorder(TreeNode<T>* p) {
  if(p == NULL) {
    cout << "NULL";
    return;
  }
  inorder(p->left);
  cout << p->element() << "/" << p->result() << "  ";
  inorder(p->right);
}

/*
 * Class Tree
 */
template<typename T>
class Tree {
  public:
    /// pointer to root node
    TreeNode<T>* root;

   /*
    * Tree constructor
    */
   Tree(double k, T f) : root(new TreeNode<T>(k, f)) { }
   /// prints tree in 'inorder' order
   void print_inorder() { inorder(root); }
};

/*
 * Class QueueNode
 * template class
 *
 * Used to store elements of priority queue
 */
template<typename T>
class QueueNode {
  public:
    /// pointer to a tree
    Tree<T>* root;
    /// pointer to next element in queue
    QueueNode<T> *next;

    /*
     * QueueNode constructor
     *
     * @param T pointer to a tree
     */
    QueueNode(Tree<T>* tree) : root(tree), next(NULL) { }
};

/*
 * Class Queue
 * template class
 *
 * Implements priority queue
 */
template<typename T>
class Queue {
  public:
    /// pointer to head of a queue
    QueueNode<T>* Head;

    /*
     * Queue constructor
     */
    Queue() : Head(new QueueNode<T>(NULL)) { }
    /// returns first element position (pointer to Head)
    QueueNode<T>* first() { return Head; }
    /// returns last element in a queue
    QueueNode<T>* last() {
      QueueNode<T> *temp = Head;
      while(temp->next != NULL) temp = temp->next;
      return temp;
    }
    /*
     * Returns a key from a tree
     *
     * @param p pointer to a queue node
     */
    double element(QueueNode<T>* p) { return p->next->root->root->key; }
    /// prints elements of a queue
    void print() {
      QueueNode<T> *temp;
      temp = Head;
      while(temp != last()) {
        cout << element(temp) << " ";
        temp = temp->next;
      }
      cout << '\n';
    }
    /*
     * Inserts a tree to a queue
     *
     * @param T pointer to a tree
     */
    void insert (Tree<T>* tree) {
      QueueNode<T>* temp = new QueueNode<T>(tree);

      QueueNode<T> *temp1;
      temp1 = Head;
      while((temp1 != last()) && (element(temp1) < temp->root->root->key)) {
        temp1 = temp1-> next;
      }
      QueueNode<T> *p;
      p = temp1 -> next;
      temp1->next = temp;
      temp1->next->next = p;
    }
    /// delete tree from the queue
    void qdelete() {
      QueueNode<T>* temp = Head->next;
      if(temp != NULL) {
        Head->next = temp->next;
        delete temp;
      }
    }
};

/*
 * Prints a matrix
 * template function
 *
 * @param M matrix to print
 * @param n number of rows to print
 * @param m number of columns to print
 */
template<typename T>
void print_matrix(Matrix<T> M, int n, int m) {
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++)
      cout << std::fixed << std::setprecision(12) << M[i][j] << "  ";
    cout << '\n';
  }
}

/*
 * Prints a vector
 * template function
 *
 * @param r row to print
 * @param n number of elements to print
 */
template<typename T>
void print_vector(Row<T> r, int n) {
  for(int i = 0; i < n; i++)
    cout << std::setprecision(12) << r[i] << " ";
}

/*
 * Merges two trees
 * template function
 *
 * @param p pointer to a tree
 * @param q pointer to a tree
 * @return pointer to new tree which is join
 */
template<typename T>
Tree<T>* merge(Tree<T>* p, Tree<T>* q) {
  Tree<T>* temp;
  float v;
  v = p->root->result() + q->root->result();
  // temp = new Tree<T>(abs(log2(abs(v))), v);
  // temp = new Tree<T>(abs(log2(p->root->element() + q->root->element())), v);
  temp = new Tree<T>(p->root->element() + q->root->element(), v);
  temp->root->left = p->root;
  temp->root->right = q->root;
  return temp;
}

void print_debug(string msg);
void print_result(string label, double error);
