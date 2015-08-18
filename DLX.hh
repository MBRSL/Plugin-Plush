/*
 * Original code from http://uaasoftware.com/blog/demos-tutorials/dancing-links-solving-sudoku/
 */
#include "FilteredTriMesh.hh"

#include <queue>

#include <QObject>

class DLX_Node {
public:
    DLX_Node* row_header_;
    DLX_Node* col_header_;
    
    DLX_Node* left_;
    DLX_Node* right_;
    DLX_Node* up_;
    DLX_Node* down_;

    enum Type {Root, Row_header, Column_header, Node} type_;
    int row_;
    int col_;
    size_t row_count_;
    
    bool covered = false;
    
    Patch_boundary *patch_;
    
    DLX_Node() :
    DLX_Node(Node, 0, 0)
    {
    }
    DLX_Node(Type type, int row, int col) :
    row_header_(nullptr),
    col_header_(nullptr),
    left_(nullptr),
    right_(nullptr),
    up_(nullptr),
    down_(nullptr),
    type_(type),
    row_(row),
    col_(col),
    row_count_(0),
    patch_(nullptr)
    {
    }
};

//template <class Derived>
//class DLX_Iterator {
//public:
//    DLX_Iterator(DLX_Node *node) :
//    current_node_(node), is_first_iteration_(true) {}
//
//    DLX_Node* operator*() const {
//        return current_node_;
//    }
//    template <class Iter>
//    bool operator==(const Iter &other) const {
//        return !is_first_iteration_ && this->operator*() == *other;
//    }
//    template <class Iter>
//    bool operator!=(const Iter &other) const {
//        return is_first_iteration_ || this->operator*() != *other;
//    }
//    Derived& operator++() {
//        return static_cast<const Derived*>(this)->operator++();
//    }
//protected:
//    DLX_Node *current_node_;
//    bool is_first_iteration_;
//};
//
//class DLX_Down_Iterator : public DLX_Iterator<DLX_Down_Iterator> {
//public:
//    DLX_Down_Iterator(DLX_Node *header) : DLX_Iterator (header->down_) {}
//    DLX_Down_Iterator& operator++() {
//        is_first_iteration_ = false;
//        current_node_ = current_node_->down_;
//        return *this;
//    }
//};
//
//class DLX_Up_Iterator : public DLX_Iterator<DLX_Up_Iterator> {
//public:
//    DLX_Up_Iterator(DLX_Node *header) : DLX_Iterator (header->up_) {}
//    DLX_Up_Iterator& operator++() {
//        is_first_iteration_ = false;
//        current_node_ = current_node_->up_;
//        return *this;
//    }
//};
//
//class DLX_Right_Iterator : public DLX_Iterator<DLX_Right_Iterator> {
//public:
//    DLX_Right_Iterator(DLX_Node *header) : DLX_Iterator (header->right_) {}
//    DLX_Right_Iterator& operator++() {
//        is_first_iteration_ = false;
//        current_node_ = current_node_->right_;
//        return *this;
//    }
//};
//
//class DLX_Left_Iterator : public DLX_Iterator<DLX_Left_Iterator> {
//public:
//    DLX_Left_Iterator(DLX_Node *header) : DLX_Iterator (header->left_) {}
//    DLX_Left_Iterator& operator++() {
//        is_first_iteration_ = false;
//        current_node_ = current_node_->left_;
//        return *this;
//    }
//};
//
//template <class Iter>
//class DLX_Iterator_Range {
//public:
//    Iter begin() {
//        return Iter(node_);
//    }
//    Iter end() {
//        return Iter(node_);
//    }
//    
//    DLX_Iterator_Range(DLX_Node *node) : node_(node) {}
//    
//private:
//    DLX_Node *node_;
//};

class DLX : public QObject
{
    Q_OBJECT
signals:
    void setJobState(int val);

private:
    size_t nRow_;
    size_t nCol_;
    
    /// maximum results recorded
    const static size_t max_results = 20;
    /// +1 if the result is not better than previous result for each calculation
    size_t skip_counter;
    /// The search terminates if skip counter reach this number
    const static size_t max_skip_num = 100000;
    
    /*
     * Row headers are not treated as parts of the DLX, they just point to the starting/ending nodes.
     * The up/down/left/right pointer of inner nodes will not point to these headers.
     */
    DLX_Node *row_headers_;
    /*
     * Unlike row headers, column headers are parts of the DLX. They are at the -1th row of each column
     * The up/down/left/right pointer of inner nodes will not point to these headers.
     */
    DLX_Node *col_headers_;
    /*
     * Can be treated as an entrance point of DLX.
     * Only point to the first/last column headers, not part of DLX.
     */
    DLX_Node root_;
    
    std::priority_queue<    std::pair< double, std::vector<DLX_Node*> >,
                            std::vector< std::pair< double, std::vector<DLX_Node*> > >,
                            std::less< std::pair< double, std::vector<DLX_Node*> > > > results_;

    std::vector<DLX_Node*> intermediate_result_;
    size_t nMerged_subMeshes_;
    bool finished_;
    
public:
    DLX(std::vector<Patch_boundary> &subset, size_t nSubMeshes);
    ~DLX();
    
    // --> DLX Algorithm functions
    DLX_Node *choose_column(void);
    void cover(DLX_Node* col_header);
    void uncover(DLX_Node* col_header);
    void print_solution(void);
    void search(size_t k);
    
    void print_covered();
    
    std::vector< std::vector<int> > get_results();
};