export function transpose(matrix) {
  return matrix.reduce(
    (prev, next) => next.map((item, i) => (prev[i] || []).concat(next[i])),
    []
  );
}
export function transpose2(matrix) {
  return matrix.reduce(
    (prev, next) => Object.values(next).map((item, i) => (prev[i] || []).concat(item)),
    []
  );
}
