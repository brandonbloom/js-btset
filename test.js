let bts = require('.');

let xs = bts.empty;

xs = xs.conj(5);
xs = xs.conj(10);
xs = xs.conj(10);
xs = xs.conj(15);
console.log(xs.count());

for (let i = 100; i < 500; i++) {
  xs = xs.conj(i);
}
console.log(xs.count());

xs = xs.disj(10);
console.log(xs.count());

for (let i = 200; i < 300; i++) {
  xs = xs.disj(i);
}
//console.log(JSON.stringify(xs, null, 2));
console.log(xs.count());


console.log(xs.lookup(123, 123));
console.log(xs.lookup(890, 'not found'));
