#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <time.h>
#include "my_util.h"
#include <cassert>
#include <vector>
#include <cmath>

template<class T>
class tree_node{
public:
    T x;
    long y;
    tree_node<T> * left = nullptr;
    tree_node<T> * right = nullptr;
    tree_node<T> * parent = nullptr;
    tree_node<T> * prev = nullptr;
    tree_node<T> * next = nullptr;
    int height = 1;

    tree_node(){

    }
};

long total_inserts = 0;

template<class T>
class ch_tree{
public:


    tree_node<T> * rootUpper = nullptr;
    tree_node<T> * rootLower = nullptr;

    ch_tree(){

    }

    ~ch_tree(){

    }
    void double_vals(){
        double_vals(rootUpper);
        double_vals(rootLower);
    }

    void double_vals(tree_node<T> * r){
        if(r!= nullptr){
            r->y = 2 * r->y;
            double_vals(r->left);
            double_vals(r->right);
        }
    }

    void rec_delete(tree_node<T> * r){
        if(r!= nullptr){
            rec_delete(r->left);
            rec_delete(r->right);
            delete r;
        }
    }

    void reset(){
        rec_delete(rootUpper);
        rec_delete(rootLower);
        rootUpper = nullptr;
        rootLower = nullptr;
    }

    
    void remove(long x){
        rootUpper = remove(x, rootUpper);
        rootLower = remove(x, rootLower);
    }
    
    void insert(long x, long y){
        total_inserts++;

        //Unused
        tree_node<T> * next_node;
        
        tree_node<T> * insertedUpper; 
        rootUpper = insert(x, y, rootUpper, &insertedUpper);
        rootUpper = fix(insertedUpper, rootUpper, true, &next_node);
        
        tree_node<T> * insertedLower;
        rootLower = insert(x, y, rootLower, &insertedLower);
        rootLower = fix(insertedLower, rootLower, false, &next_node);

        //isConvex();
    }

    void intact(){
        intact(rootUpper);
        intact(rootLower);
    }

    void isConvex(){
        isConvex(rootUpper, true);
        isConvex(rootLower, false);
    }

    void printInline(){
        std::cout<<"upperHull"<<std::endl;
        printInline(rootUpper);
        std::cout<<"lowerHull"<<std::endl;
        printInline(rootLower);
    }

    void remove_seg(long x_from, long x_to){
        return;
        assert(x_from <= x_to);
        DBGP(remove_seg)
        DBGPB(x_from, x_to);
        if (rootLower == nullptr || rootUpper == nullptr){
            return;
        }
        rootUpper = remove_seg(rootUpper, x_from, x_to);
        rootLower = remove_seg(rootLower, x_from, x_to);
    }


    //TODO this can obv be optimized
    static tree_node<T> * remove_seg(tree_node<T> * root, long x_from, long x_to){
        tree_node<T> * tn = find_geq(x_from, root);
        //DBGPA(tn->x)
        while(true){
            if (tn != nullptr){
                if(tn->x <= x_to){
                    tree_node<T> * tdel = tn;
                    tn = tn->next;
                    root = remove(tdel->x, root);
                }
                else{
                    return root;
                }
            }
            else{
                return root;
            }
        }
    }

    static void isConvex(tree_node<T> * node, bool upper){
        if (node != nullptr){
            if(node->prev != nullptr && node->next != nullptr){
                double cr = cross(node->prev->x, node->prev->y, node->x, node->y, node->next->x, node->next->y);
                if( upper ? (cr >= 0) : (cr <= 0) ){
                    DBGPB(cr, upper)
                    DBGPB(node->prev->x, node->prev->y)
                    DBGPB(node->x - node->prev->x, node->y - node->prev->y)
                    DBGPB(node->next->x - node->prev->x, node->next->y - node->prev->y)
                    assert(false);

                }
                else {
                    isConvex(node->right, upper);
                    isConvex(node->left, upper);
                }
            }
        }
    }
    /*
        returns width and updates a,b
    */
    std::pair<tree_node<T>*, tree_node<T>*> get_min_width(double * width, double * a, double * b){
        if(rootLower == nullptr){
            *width = 0;
            *a = 0;
            *b = 0;
            return std::pair<tree_node<T>*, tree_node<T>*>(nullptr, nullptr);
        }
        double w_upper;
        double w_upper_a;
        double w_upper_b;
        auto pp_upper = get_min_width(rootUpper, rootLower, true, &w_upper, &w_upper_a, &w_upper_b);

        double w_lower;
        double w_lower_a;
        double w_lower_b;
        auto pp_lower = get_min_width(rootLower, rootUpper, false, &w_lower, &w_lower_a, &w_lower_b);
        
        /*DBGPB(w_lower,w_upper)
        if(pp_lower.first != nullptr){
            DBGPB(pp_lower.first->x, pp_lower.first->y)
            DBGPB(pp_lower.second->x, pp_lower.second->y)
        }
        if(pp_upper.first != nullptr){
            DBGPB(pp_upper.first->x, pp_upper.first->y)
            DBGPB(pp_upper.second->x, pp_upper.second->y)
        }*/

        if(w_lower < w_upper){
            * a = w_lower_a;
            if (w_lower != INFINITY){
                * b = w_lower_b + (0.5 * w_lower);
                *width = w_lower;
            }
            else{
                * b = w_lower_b;
                *width = 0;
            }
            return pp_lower;
        }
        else{
            * a = w_upper_a;
            if (w_lower != INFINITY){
                * b = w_upper_b - (0.5 * w_upper);
                *width = w_upper;
            }
            else{
                * b = w_upper_b;
                *width = 0;
            }
            
            return pp_upper;
        }
    }

    //if upper == true
    //root=>upper, root2=>lower
    static std::pair<tree_node<T>*, tree_node<T>*> get_min_width(tree_node<T> * root, tree_node<T> * root2, bool upper, double * width, double * slope, double * offset){
        assert(root != nullptr);
        assert(root2 != nullptr);

        if(root->next == nullptr){
            if(root->left != nullptr){
                return get_min_width(root->left, root2, upper, width, slope, offset);
            }
            else{
                *width = INFINITY;
                *slope = 0;
                *offset = root->y;
                return std::pair<tree_node<T>*, tree_node<T>*>(root, root);
            }
        }

        bool doLeft = false;
        bool doRight = false;

        tree_node<T> * p = findAntipodalPoint(root2, (root->next->x - root->x), (root->next->y - root->y), !upper);
        double w_a;
        double w_b;
        double w = get_width(root, root->next, p, &w_a, &w_b);
        
        tree_node<T> * pleft;
        double wleft;
        if(root->left != nullptr){
            pleft = findAntipodalPoint(root2, (root->left->next->x - root->left->x), (root->left->next->y - root->left->y), !upper);
            double wleft_a;
            double wleft_b;
            wleft = get_width(root->left, root->left->next, pleft, &wleft_a, &wleft_b);
        }

        tree_node<T> * pright;
        double wright;
        if(root->right != nullptr){
            pright = findAntipodalPoint(root2, (root->right->x - root->right->prev->x), (root->right->y - root->right->prev->y), !upper);
            double wright_a;
            double wright_b;
            wright = get_width(root->right->prev, root->right, pright, &wright_a, &wright_b);
        }

        if(root->left != nullptr && root->right != nullptr){
            if((wleft <= w && w < wright) || (wleft < w && w <= wright)){
                doLeft = true;
            }
            else if((w < wleft && wright <= w) || (w <= wleft && wright < w)){
                doRight = true;                
            }
            else if(w <= wleft && w <= wright){
                double wl;
                double wl_a;
                double wl_b;
                std::pair<tree_node<T>*, tree_node<T>*> pl = get_min_width(root->left, root2, upper, &wl, &wl_a, &wl_b);

                double wr;
                double wr_a;
                double wr_b;
                std::pair<tree_node<T>*, tree_node<T>*> pr = get_min_width(root->right, root2, upper, &wr, &wr_a, &wr_b);

                if(w < wl && w < wr){
                    *width = w;
                    *slope = w_a;
                    *offset = w_b;
                    return std::pair<tree_node<T>*, tree_node<T>*>(p,root);
                }
                else if(wl < wr){
                    *width = wl;
                    *slope = wl_a;
                    *offset = wl_b;
                    return pl;
                }
                else{
                    *width = wr;
                    *slope = wr_a;
                    *offset = wr_b;
                    return pr;
                }
            }
            else{
                //most likely a floating point error occurred somewhere
                //DBGPC(wleft, w, wright)
                if(wleft < wright){
                    doLeft = true;
                }
                else if (wleft > wright){
                    doRight = true;
                }
                else{
                    //panic
                    *width = w;
                    *slope = w_a;
                    *offset = w_b;
                    return std::pair<tree_node<T>*, tree_node<T>*>(p,root);
                }

            }
        }
        else if(root->left == nullptr && root->right != nullptr){
            doRight = true;
        }
        else if(root->left != nullptr && root->right == nullptr){
            doLeft = true;
        }
        else{
            *width = w;
            *slope = w_a;
            *offset = w_b;
            return std::pair<tree_node<T>*, tree_node<T>*>(p,root);
        }

        if(doRight){
            double wtmp;
            double wtmp_a;
            double wtmp_b;
            std::pair<tree_node<T>*, tree_node<T>*> ptmp = get_min_width(root->right, root2, upper, &wtmp, &wtmp_a, &wtmp_b);
            if (wtmp < w){
                *width = wtmp;
                *slope = wtmp_a;
                *offset = wtmp_b;
                return ptmp;
            }
            else{
                *width = w;
                *slope = w_a;
                *offset = w_b;
                return std::pair<tree_node<T>*, tree_node<T>*>(p,root);
            }
        }
        else if(doLeft){
            double wtmp;
            double wtmp_a;
            double wtmp_b;
            std::pair<tree_node<T>*, tree_node<T>*> ptmp = get_min_width(root->left, root2, upper, &wtmp, &wtmp_a, &wtmp_b);
            if (wtmp < w){
                *width = wtmp;
                *slope = wtmp_a;
                *offset = wtmp_b;
                return ptmp;
            }
            else{
                *width = w;
                *slope = w_a;
                *offset = w_b;
                return std::pair<tree_node<T>*, tree_node<T>*>(p,root);
            }
        }
        assert(false);
    }

    static double get_width(tree_node<T> * fr, tree_node<T> * to, tree_node<T> * p, double * sl, double * off){
        long x_fr = fr->x;
        long y_fr = fr->y;
        long x_to = to->x;
        long y_to = to->y;
        long x = p->x;
        long y = p->y;
        
        double dx = x_to - x_fr;
        double dy = y_to - y_fr;

        *sl = dy/dx;
        *off = y_fr - (x_fr * *sl);
   
        double test_pos = y;
        double test_eval = (x * *sl) + *off;
        return std::abs(test_pos - test_eval);

        //double test_eval_x = (x * sl) + off;
        //double test_eval_y = (y - off) - sl;
        //double vx = ((y - test_eval_x)/2) * ((y - test_eval_x)/2);
        //double vy = ((x - test_eval_y)/2) * ((x - test_eval_y)/2);
        //return std::abs(std::sqrt(vx + vy));
    }

    static inline tree_node<T> * findAntipodalPoint(tree_node<T> * root, long x, long y, bool upper){
        tree_node<T> * prev_root = nullptr;
        while(root != nullptr){
            if (root->prev != nullptr){
                
                double cr = cross(root->x - root->prev->x, root->y - root->prev->y, x, y);

                if(upper ? (cr > 0) : (cr < 0)){
                    root = root->left;
                }
                else{
                    if(prev_root != nullptr ? prev_root->x < root->x : true){
                        prev_root = root;
                    }
                    root = root->right;
                }
            }
            else{
                if(root->left == nullptr){
                    prev_root = root;
                    root = root->right;
                }
                else{
                    root = root->left;
                }
            }
        }
        return prev_root;
    }

    static inline tree_node<T> * fix (tree_node<T> * inserted, tree_node<T> * root, bool upper, tree_node<T> ** next_node){
        assert(inserted != nullptr);
        
        *next_node = inserted->next;

        if(inserted->prev != nullptr && inserted->next != nullptr){
            double cr = cross(inserted->prev->x, inserted->prev->y, inserted->x, inserted->y, inserted->next->x, inserted->next->y);
            if(upper ? (cr > 0) : (cr < 0)){
                root = remove(inserted->x, root);
                return root;
            }
        }
        if (inserted->next != nullptr){
            tree_node<T> * curr = inserted;
            tree_node<T> * ncurr = curr->next;
            while(ncurr->next != nullptr){
                tree_node<T> * nncurr = ncurr->next;
                double cr = cross(curr->x, curr->y, ncurr->x, ncurr->y, nncurr->x, nncurr->y);
                if(upper ? (cr > 0) : (cr < 0)){
                    root = remove(ncurr->x, root);
                    ncurr = nncurr;
                    curr = ncurr->prev;
                }
                else{
                    break;
                }
            }
            inserted = curr;
            *next_node = inserted->next;
        }
        if (inserted->prev != nullptr){
            tree_node<T> * curr = inserted;
            tree_node<T> * ncurr = curr->prev;
            while(ncurr->prev != nullptr){
                tree_node<T> * nncurr = ncurr->prev;
                double cr = cross(nncurr->x, nncurr->y, ncurr->x, ncurr->y, curr->x, curr->y);
                if(upper ? (cr > 0) : (cr < 0)){
                    root = remove(ncurr->x, root);
                    ncurr = curr->prev;
                }
                else{
                    break;
                }
            }
        }
        return root;
    }

    static inline int get_height(tree_node<T> * tree){
        if (tree == nullptr){
            return 0;
        }
        else {
            return tree->height;
        }
    }

    static void update_height(tree_node<T> * tree){
        if (tree != nullptr){
            tree->height = std::max(get_height(tree->left) + 1, get_height(tree->right) + 1);
        }
    }

    static tree_node<T> * rotate_right(tree_node<T> * next){
        tree_node<T> * p = next->parent;
        tree_node<T> * b = next->left->right;
        tree_node<T> * n = next;
        next = next->left;
        next->parent = p;
        if (p != nullptr){
            if (p->left == n){
                p->left = next;
            }else{
                p->right = next;
            }
        }
        next->right = n;
        n->parent = next;
        n->left = b;
        if (b !=nullptr){
            b->parent = n;
        }
        update_height(next->left);
        update_height(next->right);
        return next;
    }

    static tree_node<T> * rotate_left(tree_node<T> * next){
        tree_node<T> * p = next->parent;
        tree_node<T> * b = next->right->left;
        tree_node<T> * n = next;
        next = next->right;
        next->parent = p;
        if (p != nullptr){
            if (p->left == n){
                p->left = next;
            }else{
                p->right = next;
            }
        }
        next->left = n;
        n->parent = next;
        n->right = b;
        if (b !=nullptr){
            b->parent = n;
        }
        update_height(next->left);
        update_height(next->right);
        return next;
    }

    static tree_node<T> * rebalance(tree_node<T> * next){
        if (get_height(next->left) > get_height(next->right) + 1 ){
            //rebalance the left side
            if (get_height(next->left->right) > get_height(next->left->left)){
                rotate_left(next->left);
            }
            next = rotate_right(next);
        }
        else if (get_height(next->right) > get_height(next->left) + 1 ){
            //rebalance the right side
            if (get_height(next->right->left) > get_height(next->right->right)){
                rotate_right(next->right);
            }
            next = rotate_left(next);
        }
        return next;
    }

    static tree_node<T> * insert(long x, long y, tree_node<T> * tree, tree_node<T> ** inserted){
        if (tree == nullptr){
            tree = new tree_node<T>();
            *inserted = tree;
            tree->x = x;
            tree->y = y;
            return tree;
        }
        tree_node<T> * next = tree;
        while (true){
            if (x < next->x){
                //left
                if (next->left != nullptr){
                    next = next->left;
                }
                else{
                    next->left = new tree_node<T>();
                    *inserted = next->left;
                    next->left->parent = next;
                    next = next->left;
                    //Linked list insert
                    next->next = next->parent;
                    next->prev = next->parent->prev;
                    if(next->prev != nullptr){
                        next->prev->next = next;
                    }
                    next->parent->prev = next;
                    break;
                }           
            }else if (x > next->x){
                //right
                if (next->right != nullptr){
                    next = next->right;
                }
                else{
                    next->right = new tree_node<T>();
                    *inserted = next->right;
                    next->right->parent = next;
                    next = next->right;
                    //Linked list insert
                    next->prev = next->parent;
                    next->next = next->parent->next;
                    if(next->next != nullptr){
                        next->next->prev = next;
                    }
                    next->parent->next = next;
                    break;
                }
            }
            else{
                //found
                //assert(false);
                *inserted = next;
                next->y = y;
                return tree;
            }
        }

        next->x = x;
        next->y = y;
        next->height = 1;
        next = next->parent;
        while(next != tree){
            next = rebalance(next);
            update_height(next);
            next = next->parent;
        }
        tree = rebalance(tree);
        update_height(tree);
        return tree;
    }

    static tree_node<T> * find_leq (long x,  tree_node<T> * tree){
        tree_node<T> * next = tree;
        while (true){
            if (x < next->x){
                //left
                if (next->left != nullptr){
                    next = next->left;
                }
                else{
                    return next->prev;
                }           
            }else if (x > next->x){
                //right
                if (next->right != nullptr){
                    next = next->right;
                }
                else{
                    return next;
                }
            }
            else{
                //found
                return next;
            }
        }
    }

    static tree_node<T> * find_geq (long x,  tree_node<T> * tree){
        tree_node<T> * next = tree;
        while (true){
            if (x < next->x){
                //left
                if (next->left != nullptr){
                    next = next->left;
                }
                else{
                    return next;
                }           
            }else if (x > next->x){
                //right
                if (next->right != nullptr){
                    next = next->right;
                }
                else{
                    return next->next;
                }
            }
            else{
                //found
                return next;
            }
        }
    }

    static tree_node<T> * remove(long x, tree_node<T> * tree){
        if(tree == nullptr)
        {
            return nullptr;
        }

        if(x == tree->x)
        {
            if(tree->left == nullptr || tree->right == nullptr)
            {
                //Linked list delete
                if(tree->next != nullptr){
                    tree->next->prev = tree->prev;
                }
                if(tree->prev != nullptr){
                    tree->prev->next = tree->next;
                }
                if(tree->left == nullptr && tree->right == nullptr)
                {
                    delete tree;
                    return nullptr;
                }
                else if(tree->left == nullptr)
                {
                    //Tree delete
                    tree_node<T> * tmp = tree->right;
                    tmp->parent = tree->parent;
                    delete tree;
                    tree = tmp;
                }
                else if(tree->right == nullptr)
                {
                    //Tree delete
                    tree_node<T> * tmp = tree->left;
                    tmp->parent = tree->parent;
                    delete tree;
                    tree = tmp;
                }
            }
            else
            {
                tree_node<T> *t = tree->left;
                while(t->right != nullptr)
                {
                    t = t->right;
                }
                
                long temp = tree->x;
                tree->x = t->x;
                t->x = temp;

                long tempy = tree->y;
                tree->y = t->y;
                t->y = tempy;

                t->next = tree->next;
                tree->prev = t->prev;
                tree->next = t;
                t->prev = tree;

                if(tree->prev != nullptr){
                    tree->prev->next = tree;
                }
                if(t->next != nullptr){
                    t->next->prev = t;
                }

                tree->left = remove(x, tree->left);
            }
        }
        else if(x < tree->x)
        {
            tree->left = remove(x, tree->left);
        }
        else{
            tree->right = remove(x, tree->right);
        }
        // Rebalance tree
        tree_node<T> * ret = rebalance(tree);
        update_height(ret);
        return ret;
    }

    static std::string print(tree_node<T> * tree){
        if (tree == nullptr){
            return "_";
        }
        return "(x=" + std::to_string(tree->x) + ",y=" + std::to_string(tree->y) + ",h=" + std::to_string(tree->height) + " {" + print(tree->left) + " : " + print(tree->right) + "})";
    }

    static void printInline(tree_node<T> * tree){
        if (tree == nullptr){
            std::cout << "Empty" << std::endl;
            return;
        }
        while(tree->left != nullptr){
            tree = tree->left;
        }
        while(tree != nullptr){
            std::cout << "x=" << tree->x << "y=" << tree->y << std::endl;
            tree = tree->next;
        }
    }

    static std::string printshort(tree_node<T> * tree, int tabbing){
        std::string tabs = "";
        for(int i = 0; i < tabbing; i++){
            tabs += "\t";
        }
        if (tree == nullptr){
            return tabs + "nullptr";
        }
        return tabs + "x=" + std::to_string(tree->x) + "y=" + std::to_string(tree->y) + ",h=" + std::to_string(tree->height) + ",prev=" + std::to_string(tree->prev != nullptr ? tree->prev->x : -1) + ",next=" + std::to_string(tree->next != nullptr ? tree->next->x : -1) + ",parent=" + std::to_string(tree->parent != nullptr ? tree->parent->x : -1) + "\n" + printshort(tree->left, tabbing+1) + "\n" + printshort(tree->right,tabbing+1) + "\n";
    }

    static void intact(tree_node<T> * tree){
        if (tree == nullptr){
            return;
        }
        assert(std::abs(get_height(tree->left) - get_height(tree->right)) <= 1);

        if(tree->next != nullptr){
            assert(tree->next->prev == tree);
            assert(tree->next->x > tree->x);
        }
        if(tree->prev != nullptr){
            assert(tree->prev->next == tree);
            assert(tree->prev->x < tree->x);
        }

        if(tree->left != nullptr){
            assert(tree->left->x < tree->x);
            intact(tree->left);
        }
        if(tree->right != nullptr){
            assert(tree->x < tree->right->x);
            intact(tree->right);
        }

    }

    static void assertPrint(bool assertion, tree_node<T> * tree){
        if (!assertion){
            std::cout << printshort(tree, 0) << "\n";
            assert(false);
        }
    }

    static void intactPrint(tree_node<T> * tree){
        if (tree == nullptr){
            return;
        }
        assertPrint(std::abs(get_height(tree->left) - get_height(tree->right)) <= 1, tree);

        if(tree->next != nullptr){
            assertPrint(tree->next->prev == tree, tree);
            assertPrint(tree->next->x > tree->x, tree);
        }
        if(tree->prev != nullptr){
            assertPrint(tree->prev->next == tree, tree);
            assertPrint(tree->prev->x < tree->x, tree);
        }
        if(tree->left != nullptr){
            assertPrint(tree->left->x < tree->x, tree);
            intactPrint(tree->left);
        }
        if(tree->right != nullptr){
            assertPrint(tree->x < tree->right->x, tree);
            intactPrint(tree->right);
        }
    }

    static tree_node<T>* insert_node(tree_node<T>* tree, tree_node<T>* node){
        assert(node!=nullptr);
        if (tree == nullptr){
            tree = node;
            return tree;
        }
        node->right = nullptr;
        node->left = nullptr;

        tree_node<T>* next = tree;
        while (true){
            if (node->x < next->x){
                //left
                if (next->left != nullptr){
                    next = next->left;
                }
                else{
                    next->left = node;
                    next->left->parent = next;
                    next = next->left;
                    break;
                }           
            }else if (node->x > next->x){
                //right
                if (next->right != nullptr){
                    next = next->right;
                }
                else{
                    next->right = node;
                    next->right->parent = next;
                    next = next->right;
                    break;
                }
            }
            else{
                //found
                assert(false);
                return tree;
            }
        }

        next->height = 1;
        next = next->parent;
        while(next != tree){
            next = rebalance(next);
            update_height(next);
            next = next->parent;
        }
        tree = rebalance(tree);
        update_height(tree);
        return tree;
    }

    static tree_node<T>* merge(tree_node<T>*tree1, tree_node<T>*tree2, tree_node<T>*mergeNode)
    {
        if(mergeNode == nullptr){
            return tree1;
        }

        if(tree1 == nullptr && tree2 == nullptr){
            return mergeNode;
        }
        else if(tree1 == nullptr)
        {
            tree2->parent = nullptr;

            tree_node<T>* inserted;
            tree2 = insert_node(tree2, mergeNode);
            return tree2;
        }
        else if(tree2 == nullptr)
        {
            tree1->parent = nullptr;

            tree_node<T>* inserted;
            tree1 = insert_node(tree1, mergeNode);
            return tree1;
        }
        tree1->parent = nullptr;
        tree2->parent = nullptr;

        if(get_height(tree1) - get_height(tree2) >= 2)
        {
            tree1->left = merge(tree1->left, tree2, mergeNode);
            tree1->left->parent = tree1;
            tree1 = rebalance(tree1);
            return tree1;
        }
        else if(get_height(tree2) - get_height(tree1) >= 2)
        {
            tree2->right = merge(tree1, tree2->right, mergeNode);
            tree2->right->parent = tree2;
            tree2 = rebalance(tree2);
            return tree2;
        }
        mergeNode->left = tree2;
        mergeNode->left->parent = mergeNode;
        mergeNode->right = tree1;
        mergeNode->right->parent = mergeNode;

        mergeNode = rebalance(mergeNode);
        return mergeNode;
    }

    void split(ch_tree * other, T x){
        breakTree(x, &this->rootLower, &other->rootLower);

        std::swap(rootLower, other->rootLower);
        if(this->rootLower != nullptr && other->rootLower != nullptr){
            this->rootLower->parent = nullptr;
            other->rootLower->parent = nullptr;

            tree_node<T>* tn = rootLower;
            while(tn->right != nullptr){
                tn = tn->right;
            }
            if(tn != nullptr){
                if(tn->next == nullptr){
                    tree_node<T>* tn_other = other->rootLower;
                    while(tn_other->left != nullptr){
                        tn_other = tn_other->left;
                    }
                    tn_other->prev = nullptr;
                }
                else{
                    tn->next->prev = nullptr;
                }
                tn->next = nullptr;
            }
        }

        breakTree(x, &this->rootUpper, &other->rootUpper);
        
        std::swap(rootUpper, other->rootUpper);
        if(this->rootUpper != nullptr && other->rootUpper != nullptr){
            this->rootUpper->parent = nullptr;
            other->rootUpper->parent = nullptr;

            tree_node<T>* tn = rootUpper;
            while(tn->right != nullptr){
                tn = tn->right;
            }
            if(tn != nullptr){
                if(tn->next == nullptr){
                    tree_node<T>* tn_other = other->rootUpper;
                    while(tn->left != nullptr){
                        tn_other = tn_other->left;
                    }
                    tn_other->prev = nullptr;
                }
                else{
                    tn->next->prev = nullptr;
                }
                tn->next = nullptr;
            }
        }
    }



    //Splitting tree1
    static void breakTree(T key, tree_node<T>**tree1, tree_node<T>**tree2)
    {
        if(tree1 == nullptr){
            return;
        }
        if((*tree1) == nullptr){
            return;
        }

        if(key > (*tree1)->x){
            breakTree(key, &(*tree1)->right, tree2);
            tree_node<T>*t = (*tree1)->right;
            (*tree1)->right = nullptr;

            if((*tree2) != nullptr){
                (*tree2) = merge((*tree2), (*tree1)->left, (*tree1));
            }
            else{
                (*tree2) = (*tree1)->left;
                if((*tree2) != nullptr){
                    (*tree2)->parent = nullptr;
                }

                (*tree2) = insert_node((*tree2), (*tree1));
            }
            (*tree1) = t;
        }
        else if(key < (*tree1)->x)
        {
            breakTree(key, &(*tree1)->left, tree2);
            if((*tree1)->left != nullptr && (*tree1)->right != nullptr){
                (*tree1)->right->parent = nullptr;
                (*tree1) = merge((*tree1)->right, (*tree1)->left, (*tree1));
            }
            else{
                if((*tree1)->left == nullptr)
                {
                    tree_node<T> * t = (*tree1);
                    (*tree1) = (*tree1)->right;
                    
                    if((*tree1) != nullptr){
                        (*tree1)->parent = nullptr;
                    }

                    (*tree1) = insert_node((*tree1), t);
                }
                else if((*tree1)->right == nullptr)
                {
                    tree_node<T>*t = (*tree1);
                    (*tree1) = (*tree1)->left;

                    if((*tree1) != nullptr){
                        (*tree1)->parent = nullptr;
                    }
                    
                    (*tree1) = insert_node((*tree1), t);
                }
            }
        }
        else{
            tree_node<T> *n = (*tree1);
            (*tree2) = (*tree1)->left;
            if((*tree2) != nullptr){
                (*tree2)->parent = nullptr;
            }

            (*tree1) = (*tree1)->right;
            if((*tree1) != nullptr){
                (*tree1)->parent = nullptr;
            }

            if(n->next != nullptr){
                n->next->prev = n->prev;
            }
            if(n->prev != nullptr){
                n->prev->next = n->next;
            }
            
            delete n;
        }
    }
    long getSize(){
        long ret = 0;
        ret+=getSize(rootLower);
        ret+=getSize(rootUpper);
        return ret;
    }

    long getSize(tree_node<T> * r){
        long ret = 0;
        if(r!= nullptr){
            ret = 1;
            ret += getSize(r->left);
            ret += getSize(r->right);
        }
        return ret;
    }
};

