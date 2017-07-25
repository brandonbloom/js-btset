// Based on Nikita Prokopov's btset.cljc in Datascript.

let bitTest = (x, n) => (x & (1 << n)) !== 0;

let minLen = 16;
let maxLen = 32;
let avgLen = (minLen + maxLen) >>> 1;

let levelShift;
for (let i = 31; i >= 0; i--) {
  if (bitTest(maxLen, i)) {
    levelShift = i;
    break;
  }
}

let pathMask = (1 << levelShift) - 1;
let emptyPath = 0;

let pathGet = (path, level) => {
  return pathMask & (path >>> level);
};

let pathSet = (path, level, idx) => {
  return path | (idx << level);
};

//XXX Export.
let cmp = (x, y) => (x < y ? -1 : y < x ? 1 : 0);

let binSearch = (d, arr, r, k) => {
  let l = 0;
  while (l <= r) {
    let m = (l + r) >>> 1;
    let mk = arr[m];
    if (cmp(mk, k) === d) {
      l = m + 1;
    } else {
      r = m - 1;
    }
  }
  return l;
};

let binSearchL = (arr, r, k) => binSearch(-1, arr, r, k);
let binSearchR = (arr, r, k) => binSearch(+1, arr, r, k);

let lookupExact = (arr, k) => {
  let n = arr.length;
  let idx = binSearchL(arr, n - 1, k);
  if (idx < n && cmp(arr[idx], k) === 0) {
    return idx;
  }
  return -1;
};

let lookupRange = (arr, k) => {
  let n = arr.length;
  let idx = binSearchL(arr, n - 1, k);
  if (idx === n) {
    return -1;
  }
  return idx;
};

// Array operations.

let alast = arr => arr[arr.length - 1];

let acopy = (src, srcStart, srcEnd, dst, dstStart) => {
  let n = srcEnd - srcStart;
  for (let i = 0; i < n; i++) {
    dst[i + dstStart] = src[srcStart + i];
  }
};

let cutNSplice = (arr, cutFrom, cutTo, spliceFrom, spliceTo, xs) => {
  let xsL = xs.length;
  let l1 = spliceFrom - cutFrom;
  let l2 = cutTo - spliceTo;
  let l1xs = l1 + xsL;
  let newArr = new Array(l1 + xsL + l2);
  acopy(arr, cutFrom, spliceFrom, newArr, 0);
  acopy(xs, 0, xsL, newArr, l1);
  acopy(arr, spliceTo, cutTo, newArr, l1xs);
  return newArr;
};

let splice = (arr, from, to, xs) =>
  cutNSplice(arr, 0, arr.length, from, to, xs);

let insert = (arr, idx, xs) =>
  cutNSplice(arr, 0, arr.length, idx, idx, xs);

let mergeNSplit = (a1, a2) => {
  let a1L = a1.length;
  let a2L = a2.length;
  let totalL = a1L + a2L;
  let r1L = totalL >>> 1;
  let r2L = totalL - r1L;
  let r1 = new Array(r1L);
  let r2 = new Array(r2L);
  if (a1L < r1L) {
    acopy(a1, 0, a1L, r1, 0);
    acopy(a2, 0, r1L - a1L, r1, a1L);
    acopy(a2, r1L - a1L, a2L, r2, 0);
  } else {
    acopy(a1, 0, r1L, r1, 0);
    acopy(a1, r1L, a1L, r2, 0);
    acopy(a2, 0, a2L, r2, a1L - r1L);
  };
  return [r1, r2];
};

let aequal = (a1, a1From, a1To, a2, a2From, a2To) => {
  let len = a1From - a1To;
  if (len !== a2To - a2From) {
    return false;
  }
  for (let i = 0; i < len; i++) {
    if (0 !== cmp(a1[i + a1From], a2[i + a2From])) {
      return false;
    }
  }
  return true;
};

let checkNSplice = (arr, from, to, newArr) => {
  if (aequal(arr, from, to, newArr, 0, newArr.length)) {
    return arr;
  }
  return splice(arr, from, to, newArr);
};

let amapInPlace = (f, arr) => {
  let n = arr.length;
  for (let i = 0; i < n; i++) {
    arr[i] = f(arr[i]);
  }
  return arr;
};

// Splits `arr` into arrays of size between min-len and max-len,
// trying to stick to (min+max)/2.
let arrPartitionApprox = (minLen, maxLen, arr) => {
  let chunkLen = avgLen;
  let len = arr.length;
  let acc = [];
  let pos = 0;
  while (len > 0) {
    let rest = len - pos;
    if (rest <= maxLen) {
      acc.push(arr.slice(pos));
      break;
    } else if (rest >= chunkLen + minLen) {
      acc.push(arr.slice(pos, pos + chunkLen));
      pos += chunkLen;
    } else {
      let pieceLen = rest >>> 1;
      acc.push(arr.slice(pos, pos + pieceLen));
      pos += pieceLen;
    }
  }
  return acc;
};

let isSortedDistinct = (arr) => {
  let al = arr.length;
  if (al <= 1) {
    return true;
  }
  let p = arr[0];
  for (let i = 1; i < al; i++) {
    let e = arr[i];
    if (cmp(e, p) === 0) {
      return false;
    }
  }
  return true;
};

let filterSortedDups = (arr) => {
  if (isSortedDistinct(arr)) {
    return arr;
  }
  let al = arr.length;
  let acc = [];
  let p = arr[0];
  for (let i = 0; i < al; i++) {
    let e = arr[i];
    if (0 !== cmp(e, p)) {
      acc.push(e);
    }
    p = e;
  }
  return acc;
};

let returnArr = (...arr) => arr.filter(x => x !== null);

// Nodes.

let rotate = (node, root, left, right) => {
  // Root never merges.
  if (root) {
    return returnArr(node);
  }
  // Enough keys, nothing to merge.
  if (node.length() > minLen) {
    return returnArr(left, node, right);
  }
  // Left and this can be merged to one.
  if (left && left.length() <= minLen) {
    return returnArr(left.merge(node), right);
  }
  // Right and this can be merged to one.
  if (right && right.length() <= minLen) {
    return returnArr(left, node.merge(right));
  }
  // Left has fewer nodes, redestribute with it.
  if (left && (right === null || left.length() < right.length())) {
    let nodes = left.mergeNSplit(node);
    return returnArr(nodes[0], nodes[1], right);
  }
  // Right has fewer nodes, redestribute with it.
  let nodes = node.mergeNSplit(right);
  return returnArr(left, nodes[0], nodes[1]);
};

class Node {
  constructor(keys, ptrs) {
    this.keys = keys;
    this.ptrs = ptrs;
  }

  limKey() {
    return alast(this.keys);
  }

  length() {
    return this.keys.length;
  }

  merge(next) {
    return new Node(this.keys.concat(next.keys),
                    this.ptrs.concat(next.ptrs));
  }

  mergeNSplit(next) {
    let ks = mergeNSplit(this.keys, next.keys);
    let ps = mergeNSplit(this.ptrs, next.ptrs);
    return returnArr(new Node(ks[0], ps[0]), new Node(ks[1], ps[1]));
  }

  lookup(key, notFound) {
    let idx = lookupRange(this.keys, key);
    if (idx === -1) {
      return notFound;
    }
    return this.ptrs[idx].lookup(key, notFound);
  }

  conj(key) {
    let idx = binSearchL(this.keys, this.keys.length - 2, key);
    let nodes = this.ptrs[idx].conj(key);
    if (!nodes) {
      return null;
    }
    let lims = nodes.map(n => n.limKey());
    let newKeys = checkNSplice(this.keys, idx, idx + 1, lims);
    let newPtrs = splice(this.ptrs, idx, idx + 1, nodes);
    if (newPtrs.length <= maxLen) {
      // OK as is.
      return [new Node(newKeys, newPtrs)];
    }
    // Gotta split it up.
    let middle = newPtrs.length >>> 1;
    return [
      new Node(newKeys.splice(0, middle), newPtrs.splice(0, middle)),
      new Node(newKeys.splice(middle), newPtrs.splice(middle)),
    ];
  }

  disj(key, root, left, right) {
    let idx = lookupRange(this.keys, key);
    if (idx === -1) {
      // Short-circuit, key not here.
      return null
    }
    let child = this.ptrs[idx];
    let leftChild = (idx - 1 >= 0 ? this.ptrs[idx - 1] : null);
    let rightChild = (idx + 1 < this.ptrs.length ? this.ptrs[idx + 1] : null);
    let disjned = child.disj(key, false, leftChild, rightChild);
    if (disjned === null) {
      // Short-circuit, key not here.
      return null;
    }
    let leftIdx = (leftChild ? idx - 1 : idx);
    let rightIdx = (rightChild ? idx + 2 : idx + 1);
    let lims = disjned.map(n => n.limKey());
    let newKeys = checkNSplice(this.keys, leftIdx, rightIdx, lims);
    let newPtrs = splice(this.ptrs, leftIdx, rightIdx, disjned);
    return rotate(new Node(newKeys, newPtrs), root, left, right);
  }
}

class Leaf {
  constructor(keys) {
    this.keys = keys;
  }

  limKey() {
    return alast(this.keys);
  }

  length() {
    return this.keys.length;
  }

  merge(next) {
    return new Leaf(this.keys.concat(next.keys));
  }

  mergeNSplit(next) {
    let ks = mergeNSplit(this.keys, next.keys);
    return [new Leaf(ks[0]), new Leaf(ks[1])];
  }

  lookup(key, notFound) {
    let idx = lookupExact(this.keys, key);
    return (idx == -1 ? notFound : this.keys[idx]);
  }

  conj(key) {
    let idx = binSearchL(this.keys, this.keys.length - 1, key);
    let keysL = this.keys.length;
    // Element already here?
    if (idx < keysL && 0 === cmp(key, this.keys[idx])) {
      return null;
    }
    // Splitting.
    if (keysL === maxLen) {
      let middle = (keysL + 1) >>> 1;
      if (idx > middle) {
        // New key goes to the second half.
        return [
          new Leaf(this.keys.slice(0, middle)),
          new Leaf(cutNSplice(this.keys, middle, keysL, idx, idx, [key])),
        ];
      } else {
        // New key goes to the first half.
        return [
          new Leaf(cutNSplice(this.keys, 0, middle, idx, idx, [key])),
          new Leaf(this.keys.slice(middle, keysL)),
        ];
      }
    }
    // OK as is.
    return [new Leaf(splice(this.keys, idx, idx, [key]))];
  }

  disj(key, root, left, right) {
    let idx = lookupExact(this.keys, key);
    if (idx === -1) {
      return null;
    }
    let newKeys = splice(this.keys, idx, idx + 1, []);
    return rotate(new Leaf(newKeys), root, left, right);
  }
}

// BTSet.

class BTSet {
  constructor(root, shift, cnt) {
    this._root = root;
    this._shift = shift;
    this._cnt = cnt;
  }

  equal(other) {
    if (this._cnt !== other._cnt) {
      return false;
    }
    let it = this.iterator();
    while (it.hasNext()) {
      let x = it.next();
      if (!other.contains(x)) {
        return false;
      }
    }
    return true;
  }

  conj(key) {
    let roots = this._root.conj(key);
    if (roots === null) {
      // Tree not changed.
      return this;
    }
    if (roots.length === 1) {
      // Keeping single root.
      return new BTSet(roots[0], this._shift, this._cnt + 1);
    }
    // Introducing new root.
    let lims = roots.map(n => n.limKey());
    let node = new Node(lims, roots);
    return new BTSet(node, this._shift + levelShift, this._cnt + 1);
  }

  disj(key) {
    let newRoots = this._root.disj(key, true, null, null);
    if (newRoots === null) {
      // Nothing changed, key wasn't in the set.
      return this;
    }
    let newRoot = newRoots[0];
    if (newRoot instanceof Node && newRoot.length() === 1) {
      // Root has one child, make it new root.
      let shift = this._shift - levelShift;
      return new BTSet(newRoot.ptrs[0], shift, this._cnt - 1);
    }
    // Keeping root level.
    return new BTSet(newRoot, this._shift, this._cnt - 1);
  }

  lookup(k, notFound) {
    return this._root.lookup(k, notFound);
  }

  count() {
    return this._cnt;
  }

  _keysFor(path) {
    let level = this._shift;
    let node = this._root;
    while (level > 0) {
      level -= levelShift;
      node = node.ptrs[pathGet(path, level)];
    }
    return node.keys;
  }
}

/*

;; iteration

(defn -next-path ^long [node ^long path ^long level]
  (let [idx (path-get path level)]
    (if (pos? level)
      ;; inner node
      (let [sub-path (-next-path (da/aget (.-pointers ^Node node) idx) path (- level level-shift))]
        (if (== -1 sub-path)
          ;; nested node overflow
          (if (< (inc idx) (da/alength (.-pointers ^Node node)))
            ;; advance current node idx, reset subsequent indexes
            (path-set empty-path level (inc idx))
            ;; current node overflow
            -1)
          ;; keep current idx
          (path-set sub-path level idx)))
      ;; leaf
      (if (< (inc idx) (da/alength (.-keys ^Leaf node)))
        ;; advance leaf idx
        (path-set empty-path 0 (inc idx))
        ;; leaf overflow
        -1))))

(defn next-path
  "Returns path representing next item after `path` in natural traversal order,
   or -1 if end of tree has been reached"
  ^long [^BTSet set ^long path]
  (-next-path (.-root set) path (.-shift set)))

(defn -rpath
  "Returns rightmost path possible starting from node and going deeper"
  ^long [node ^long level]
  (loop [node  node
         path  empty-path
         level level]
    (if (pos? level)
      ;; inner node
      (recur (alast (.-pointers ^Node node))
             (path-set path level (dec (da/alength (.-pointers ^Node node))))
             (- level level-shift))
      ;; leaf
      (path-set path 0 (dec (da/alength (.-keys ^Leaf node)))))))

(defn -prev-path ^long [node ^long path ^long level]
  (let [idx (path-get path level)]
    (if (pos? level)
      ;; inner node
      (let [sub-level (- level level-shift)
            sub-path  (-prev-path (da/aget (.-pointers ^Node node) idx) path sub-level)]
        (if (== -1 sub-path)
          ;; nested node overflow
          (if (>= (dec idx) 0)
            ;; advance current node idx, reset subsequent indexes
            (let [idx      (dec idx)
                  sub-path (-rpath (da/aget (.-pointers ^Node node) idx) sub-level)]
              (path-set sub-path level idx))
            ;; current node overflow
            -1)
          ;; keep current idx
          (path-set sub-path level idx)))
      ;; leaf
      (if (>= (dec idx) 0)
        ;; advance leaf idx
        (path-set empty-path 0 (dec idx))
        ;; leaf overflow
        -1))))

(defn prev-path
  "Returns path representing previous item before `path` in natural traversal order,
   or -1 if `path` was already beginning of a tree"
  ^long [^BTSet set ^long path]
  (-prev-path (.-root set) path (.-shift set)))

(defn btset-iter
  "Iterator that represents whole set"
  [^BTSet set]
  (when (pos? (node-len (.-root set)))
    (let [left   empty-path
          right  (inc (-rpath (.-root set) (.-shift set)))]
      (iter set left right))))

(deftype Iter [set ^long left ^long right keys ^long idx]
  #?@(:cljs [
    ISeqable
    (-seq [this] (when keys this))

    ISeq
    (-first [this] (iter-first this))
    (-rest [this]  (or (iter-next this) ()))

    INext
    (-next [this] (iter-next this))

    IChunkedSeq
    (-chunked-first [this] (iter-chunk this))
    (-chunked-rest  [this] (or (-chunked-next this) ()))

    IChunkedNext
    (-chunked-next  [this] (iter-chunked-next this))
             
    IReduce
    (-reduce [this f] (iter-reduce this f))
    (-reduce [this f start] (iter-reduce this f start))

    IReversible
    (-rseq [this] (iter-rseq this))]
  ))

(defn iter [set ^long left ^long right]
  (Iter. set left right (keys-for set left) (path-get left 0)))

(defn iter-first [^Iter iter]
  (when (.-keys iter)
    (da/aget (.-keys iter) (.-idx iter))))

(defn iter-next [^Iter iter]
  (let [set   (.-set iter)
        left  (.-left iter)
        right (.-right iter)
        keys  (.-keys iter)
        idx   (.-idx iter)]
    (when keys
      (if (< (inc idx) (da/alength keys))
        ;; can use cached array to move forward
        (when (< (inc left) right)
          (Iter. set (inc left) right keys (inc idx)))
        (let [left (next-path set left)]
          (when (and (not= -1 left) (< left right))
            (datascript.btset/iter set left right)))))))

(defn iter-chunk [^Iter iter]
  (let [left  (.-left iter)
        right (.-right iter)
        keys  (.-keys iter)
        idx   (.-idx iter)
        end-idx (if (= (bit-or left path-mask)
                       (bit-or right path-mask))
                  (bit-and right path-mask)
                  (da/alength keys))]
      (#?(:cljs array-chunk) keys idx end-idx)))

(defn iter-chunked-next [^Iter iter]
  (let [set   (.-set iter)
        left  (.-left iter)
        right (.-right iter)
        keys  (.-keys iter)
        idx   (.-idx iter)]
    (let [left (next-path set (+ left (- (da/alength keys) idx 1)))]
      (when (and (not= -1 left) (< left right))
        (datascript.btset/iter set left right)))))

(defn iter-rseq [^Iter iter]
  (let [set   (.-set iter)
        left  (.-left iter)
        right (.-right iter)]
    (when (.-keys iter)
      (riter set (prev-path set left) (prev-path set right)))))

(defn iter-reduce
  ([^Iter iter f]
    (if (nil? (.-keys iter))
      (f)
      (let [first (iter-first iter)]
        (if-let [next (iter-next iter)]
          (iter-reduce next f first)
          first))))
        
  ([^Iter iter f start]
    (let [set   (.-set iter)
          right (.-right iter)]
      (loop [left (.-left iter)
             keys (.-keys iter)
             idx  (.-idx iter)
             acc  start]
        (if (nil? keys)
          acc
          (let [new-acc (f acc (da/aget keys idx))]
            (cond
              (reduced? new-acc)
                @new-acc
              (< (inc idx) (da/alength keys)) ;; can use cached array to move forward
                (if (< (inc left) right)
                  (recur (inc left) keys (inc idx) new-acc)
                  new-acc)
              :else
                (let [new-left (next-path set left)]
                  (if (and (not== -1 new-left) (< new-left right))
                    (recur new-left (keys-for set new-left) (path-get new-left 0) new-acc)
                    new-acc)))))))))

;; reverse iteration

(deftype ReverseIter [set ^long left ^long right keys ^long idx]
  #?@(:cljs [
    ISeqable
    (-seq [this] (when keys this))

    ISeq
    (-first [this] (riter-first this))
    (-rest [this]  (or (riter-next this) ()))

    INext
    (-next [this] (riter-next this))

    IReversible
    (-rseq [this] (riter-rseq this))]
  ))

(defn riter [set ^long left ^long right]
  (ReverseIter. set left right (keys-for set right) (path-get right 0)))

(defn riter-first [^ReverseIter riter]
  (when (.-keys riter)
    (da/aget (.-keys riter) (.-idx riter))))

(defn riter-next [^ReverseIter ri]
  (let [set   (.-set   ri)
        left  (.-left  ri)
        right (.-right ri)
        keys  (.-keys  ri)
        idx   (.-idx   ri)]
    (when keys
      (if (>= (dec idx) 0)
        ;; can use cached array to advance
        (when (> (dec right) left)
          (ReverseIter. set left (dec right) keys (dec idx)))
        (let [right (prev-path set right)]
          (when (and (not= -1 right) (> right left))
            (riter set left right)))))))

(defn riter-rseq [^ReverseIter riter]
  (let [set   (.-set   riter)
        left  (.-left  riter)
        right (.-right riter)]
    (when keys
      (let [new-left  (if (== left -1) 0 (next-path set left))
            new-right (next-path set right)
            new-right (if (== new-right -1) (inc right) new-right)]
        (iter set new-left new-right)))))


;; distance

(defn -distance [node ^long left ^long right ^long level]
  (let [idx-l (path-get left level)
        idx-r (path-get right level)]
    (if (pos? level)
      ;; inner node
      (if (== idx-l idx-r)
        (-distance (da/aget (.-pointers ^Node node) idx-l) left right (- level level-shift))
        (loop [level level
               res   (- idx-r idx-l)]
          (if (== 0 level)
            res
            (recur (- level level-shift) (* res avg-len)))))
      (- idx-r idx-l))))

(defn distance [^BTSet set ^long path-l ^long path-r]
  (cond
    (== path-l path-r) 0
    (== (inc path-l) path-r) 1
    (== (next-path set path-l) path-r) 1
    :else (-distance (.-root set) path-l path-r (.-shift set))))

(defn est-count [^Iter iter]
  (distance (.-set iter) (.-left iter) (.-right iter)))


;; Slicing

(defn -seek
  "Returns path to first element >= key,
   or -1 if all elements in a set < key"
  [^BTSet set key]
  (loop [node  (.-root set)
         path  empty-path
         level (.-shift set)]
    (let [keys-l (node-len node)]
      (if (== 0 level)
        (let [keys (.-keys ^Leaf node)
              idx  (binary-search-l (.-comparator set) keys (dec keys-l) key)]
          (if (== keys-l idx) -1 (path-set path 0 idx)))
        (let [keys (.-keys ^Node node)
              idx  (binary-search-l (.-comparator set) keys (- keys-l 2) key)]
          (recur (da/aget (.-pointers ^Node node) idx)
                 (path-set path level idx)
                 (- level level-shift)))))))

(defn -rseek
  "Returns path to the first element that is > key.
   If all elements in a set are <= key, returns `(-rpath set) + 1`.
   It’s a virtual path that is bigger than any path in a tree"
  [^BTSet set key]
  (loop [node  (.-root set)
         path  empty-path
         level (.-shift set)]
    (let [keys-l (node-len node)]
      (if (== 0 level)
        (let [keys (.-keys ^Leaf node)
              idx  (binary-search-r (.-comparator set) keys (dec keys-l) key)]
          (path-set path 0 idx))
        (let [keys (.-keys ^Node node)
              idx  (binary-search-r (.-comparator set) keys (- keys-l 2) key)]
          (recur (da/aget (.-pointers ^Node node) idx)
                 (path-set path level idx)
                 (- level level-shift)))))))

(defn -slice [set key-from key-to]
  (let [path (-seek set key-from)]
    (when-not (neg? path)
      (let [till-path (-rseek set key-to)]
        (when (> till-path path)
          (Iter. set path till-path (keys-for set path) (path-get path 0)))))))

(defn slice
  "When called with single key, returns iterator over set that contains all elements equal to the key.
   When called with two keys (range), returns iterator for all X where key-from <= X <= key-to"
  ([set key] (slice set key key))
  ([^BTSet set key-from key-to]
    (-slice set key-from key-to)))


*/

exports.empty = new BTSet(new Leaf([]), 0, 0);
