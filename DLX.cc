/*
 * Original code from http://uaasoftware.com/blog/demos-tutorials/dancing-links-solving-sudoku/
 * Modified by MBRSL
 */
#include <stdio.h>
#include "DLX.hh"

DLX::DLX(std::vector<Patch_boundary> &subset, size_t nSubMeshes)
{
    nRow_ = subset.size();
    nCol_ = nSubMeshes;
    
    skip_counter = 0;
    finished_ = false;
    nMerged_subMeshes_ = 0;
    
    row_headers_ = new DLX_Node[nRow_];
    col_headers_ = new DLX_Node[nCol_];

    // The top most node for each column
    DLX_Node* *top_most_nodes = new DLX_Node*[nCol_];
    // Previous up_ node for each column
    DLX_Node* *up_nodes = new DLX_Node*[nCol_];
    std::fill_n(up_nodes, nCol_, nullptr);
    std::fill_n(top_most_nodes, nCol_, nullptr);
    
    for (size_t j = 0; j < nRow_; j++) {
        // We start from the first one, which contains larger patch
        Patch_boundary &patch = subset[j];
        size_t count = patch.merged_subMesh_idx.count();

        // The left_ most node in this row
        DLX_Node* left_most_nodes = nullptr;
        // The previous left_ node in this row
        DLX_Node* left_node = nullptr;
        
        for (size_t i = patch.merged_subMesh_idx.find_first();
             i < nCol_ && i != patch.merged_subMesh_idx.npos;
             i = patch.merged_subMesh_idx.find_next(i)) {
            // create a new node
            DLX_Node* node = new DLX_Node(DLX_Node::Node, j, i);
            node->row_header_ = &row_headers_[j];
            node->col_header_ = &col_headers_[i];
            node->row_count_ = count;
            
            // left_ & right_
            if (left_node != nullptr) {
                node->left_ = left_node;
                left_node->right_ = node;
            } else {
                left_most_nodes = node;
            }
            left_node = node;
            
            // up_ & down_
            if (up_nodes[i] != nullptr) {
                node->up_ = up_nodes[i];
                up_nodes[i]->down_ = node;
            } else {
                top_most_nodes[i] = node;
            }
            up_nodes[i] = node;
        }
        
        // Connect first & last column in this row
        assert(left_node != nullptr);
        left_most_nodes->left_ = left_node;
        left_node->right_ = left_most_nodes;
        
        row_headers_[j].up_ = &row_headers_[(j+nRow_-1)%nRow_];
        row_headers_[j].down_ = &row_headers_[(j+1)%nRow_];
        row_headers_[j].left_ = left_node;
        row_headers_[j].right_ = left_most_nodes;
        
        row_headers_[j].type_ = DLX_Node::Row_header;
        row_headers_[j].row_ = j;
        row_headers_[j].col_ = -1;
        row_headers_[j].row_count_ = count;
        
        row_headers_[j].patch_ = &patch;
    }
    
    // For each column, connect first & last row
    DLX_Node* prev_left_node = nullptr;
    for (size_t i = 0; i < nCol_; i++) {
        // No column can be empty, otherwise no solution
        assert(top_most_nodes[i] != nullptr && "No possible solution for this case.");
        top_most_nodes[i]->up_ = &col_headers_[i];
        up_nodes[i]->down_ = &col_headers_[i];
        
        col_headers_[i].up_ = up_nodes[i];
        col_headers_[i].down_ = top_most_nodes[i];
        col_headers_[i].left_ = &col_headers_[(i+nCol_-1)%nCol_];
        col_headers_[i].right_ = &col_headers_[(i+1)%nCol_];

        col_headers_[i].type_ = DLX_Node::Column_header;
        col_headers_[i].row_ = -1;
        col_headers_[i].col_ = i;
        col_headers_[i].row_count_ = 0;
    }
    
    delete[] up_nodes;
    delete[] top_most_nodes;
    
    //Insert root
    root_.type_ = DLX_Node::Root;
    root_.left_ = &col_headers_[nCol_-1];
    root_.right_ = &col_headers_[0];
    root_.row_ = nRow_-1;
    root_.col_ = nCol_;
    col_headers_[0].left_ = &root_;
    col_headers_[nCol_-1].right_ = &root_;
}

DLX::~DLX()
{
    for (size_t i = 0; i < nCol_; i++) {
        for (DLX_Node* row_node = col_headers_[i].down_; row_node->down_ != col_headers_[i].down_; row_node = row_node->down_) {
            delete row_node;
        }
    }
    delete[] row_headers_;
    delete[] col_headers_;
}

DLX_Node* DLX::choose_column(void) {
    return root_.right_;
}

void DLX::cover(DLX_Node* col_header) {
    col_header->right_->left_ = col_header->left_;
    col_header->left_->right_ = col_header->right_;
    
    for (DLX_Node* row_node = col_header->down_; row_node != col_header; row_node = row_node->down_) {
        row_node->covered = true;

        for(DLX_Node* right_node = row_node->right_; right_node!=row_node; right_node = right_node->right_) {
            right_node->covered = true;
            right_node->up_->down_ = right_node->down_;
            right_node->down_->up_ = right_node->up_;
        }
    }
//    print_covered();
}

void DLX::uncover(DLX_Node* col_header) {
    for (DLX_Node* row_node = col_header->up_; row_node != col_header; row_node = row_node->up_) {
        row_node->covered = false;

        for(DLX_Node* left_node = row_node->left_; left_node!=row_node; left_node = left_node->left_) {
            left_node->covered = false;
            left_node->up_->down_ = left_node;
            left_node->down_->up_ = left_node;
        }
    }
    col_header->right_->left_ = col_header;
    col_header->left_->right_ = col_header;
//    print_covered();
}

void DLX::search(size_t k) {
    typedef std::pair<std::vector<HalfedgeHandle>, double> Segment_pair;
    auto comparator = [](Segment_pair p1, Segment_pair p2) ->bool {
        return p1.second < p2.second;
    };

    emit setJobState((double)nMerged_subMeshes_/nCol_*100);
    
    if(root_.left_ == &root_ && root_.right_ == &root_) {
        double score = 0;
        for (DLX_Node *node : intermediate_result_) {
            Patch_boundary *patch = node->patch_;
            score += patch->seam_score;
        }
        
        double max_score = !results_.empty() ? results_.top().first : -1;
        if (results_.size() >= max_results && score > max_score) {
            skip_counter++;
            if (skip_counter > max_skip_num) {
                finished_ = true;
            }
        } else {
            //Valid solution!
            if (results_.size() >= max_results) {
                results_.pop();
            }
            results_.push(std::make_pair(score, intermediate_result_));
            
            skip_counter = 0;
            
            print_solution();
        }
        emit setJobState((double)skip_counter/max_skip_num*100);

        return;
    }
    if (k > nCol_) {
        printf("----------- SOMETHING WRONG -----------\n");
        finished_ = true;
        return;
    }
    DLX_Node* col_header = choose_column();
    assert(col_header->type_ != DLX_Node::Root);
    cover(col_header);
    
    for (DLX_Node* row_node = col_header->down_; !finished_ && row_node != col_header; row_node = row_node->down_) {
        // Try this row node on!
        intermediate_result_.push_back(row_node->row_header_);
        nMerged_subMeshes_ += row_node->row_count_;
        for(DLX_Node* right_node = row_node->right_; right_node!=row_node; right_node = right_node->right_) {
            cover(right_node->col_header_);
        }
        search(k+1);
        // Ok, that node didn't quite work
        for(DLX_Node* left_node = row_node->left_; left_node!=row_node; left_node = left_node->left_) {
            uncover(left_node->col_header_);
        }
        intermediate_result_.pop_back();
        nMerged_subMeshes_ -= row_node->row_count_;
    }
    
    uncover(col_header);
}

void DLX::print_covered() {
    for (size_t j = 0; j < nRow_; j++) {
        bool first_iteration = true;
        int last_col = 0;
        for(DLX_Node* right_node = row_headers_[j].right_;
            right_node != nullptr && (first_iteration || right_node!=row_headers_[j].right_);
            right_node = right_node->right_, first_iteration = false) {
            
            if (!right_node->covered) {
                for (int a = last_col; a < right_node->col_; a++) {
                    printf("   ");
                }
                printf("%3d", right_node->col_);
                last_col = right_node->col_+1;
            }
        }
        if (last_col != 0) {
            printf("\n");
        }
    }
    printf("-------------------------\n");
}

void DLX::print_solution(void) {
#ifndef NDEBUG
    printf("----------- SOLUTION FOUND -----------\n");
    for(DLX_Node *row_header : intermediate_result_) {
        printf("Row %d:", row_header->row_);
        bool first_iteration = true;
        for (DLX_Node* right_node = row_header->right_;
             right_node != nullptr && (first_iteration || right_node != row_header->right_);
             right_node = right_node->right_, first_iteration = false) {
            printf("%d ", right_node->col_);
        }
        printf("\n");
    }
    printf("---------------------\n");
#endif
}

std::vector< std::vector<int> > DLX::get_results() {
    std::vector< std::vector<int> > results;
    while (!results_.empty()) {
        auto &pair = results_.top();
        
        std::vector<int> list;
        for (DLX_Node* node : pair.second) {
            list.push_back(node->row_);
        }
        results.push_back(list);
        
        results_.pop();
    }
    return results;
}
