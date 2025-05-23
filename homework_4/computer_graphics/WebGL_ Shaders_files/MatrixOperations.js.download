/*
  Matrix operations for computer graphics.
  (c) 2023-2025 Eugen Croitoru.
  Licensed under the GPLv3+.
*/

(function(globalName, root, factory ) {
    //Asynchronous Module Definition
    if (typeof define === 'function' && define.amd) { 
      define( [], factory );
    }
    //NodeJS
    else if ( typeof exports === 'object' ) {
      module.exports = factory();
    }
    //browser
    else{
      root[globalName] = factory();
    }
  }('MatrixOperations', typeof window !== "undefined" ? window : undefined, function(){
    //this constant is supposed to not be changed.
    const matrixSize = 4;
  
    function newMatrix(ret) {
      ret = ret || new Float32Array(matrixSize * matrixSize);
      return ret;
    }
  
    function newVector(ret) {
      ret = ret || new Float32Array(matrixSize);
      return ret;
    }
    
    function newIdentityMatrix(ret) {
      ret = ret || newMatrix();
      for(var p = 0; p < matrixSize; ++p)
        ret[p + p * matrixSize] = 1;
      return ret;
    }
    
    function add(a, b, ret) {
      ret = ret || newMatrix();
      var pos = 0;
      //column-major
      for(var r = 0; r < matrixSize; ++r)
        for(var c = 0; c < matrixSize; ++c) {
      ret[pos] = a[pos] + b[pos];
      ++pos;
        }
      return ret;
    }
  
    function mulMatrices(a, b, ret) {
      ret = ret || newMatrix();
      //column-major
      var pos = 0;
      for(var r = 0; r < matrixSize; ++r)
        for(var c = 0; c < matrixSize; ++c) {
      pos = c * matrixSize + r;
      ret[pos] = 0;
      for(var t = 0; t < matrixSize; ++t) {
        ret[pos] += a[t * matrixSize + r] * b[c * matrixSize + t];
      }
        }
      return ret;
    }
  
    function mulMatrixVector(a, v, ret) {
      ret = ret || newVector();
      for(var r = 0; r < matrixSize; ++r) {
        ret[r] = 0;
        for(var t = 0; t < matrixSize; ++t) {
      ret[r] += a[t * matrixSize + r] * v[t];
        }
      }
    }
  
    function mulVectorMatrix(v, b, ret) {
      ret = ret || newVector();
      var pos = 0;
      for(var c = 0; c < matrixSize; ++c) {
        pos = c * matrixSize;
        ret[pos] = 0;
        for(var t = 0; t < matrixSize; ++t) {
      ret[pos] += a[t * matrixSize] * b[c * matrixSize + t];
        }
      }
      return ret;
    }
  
    function mulVectorsInner(v0, v1, ret) {
      ret = ret || Float32;
      for(var pos = 0; pos < matrixSize; ++c)
        ret += v0[pos] * v1[pos];
      return ret;
    }
    
    function mul(a, b, ret) {
      if(a.length == matrixSize * matrixSize && b.length == matrixSize * matrixSize)
        return mulMatrices(a, b);
      if(a.length == matrixSize * matrixSize && b.length == matrixSize)
        return mulMatrixVector(a, b);
      if(a.length == matrixSize && b.length == matrixSize * matrixSize)
        return mulVectorMatrix(a, b);
      if(a.length == matrixSize && b.length == matrixSize)
        return mulVectorsInner(a, b);
    }
  
    function projectionMatrixFrustum(left, right, bottom, top, near, far, ret) {
      ret = ret || newMatrix();
      var lx = right - left; 
      var ly = top - bottom;
      var lz = far - near;
          if(lx === 0 || ly === 0 || lz === 0)
              return newIdentityMatrix();
          if(near <= 0 || far <= 0) {
              console.log("Error: Projection Matrix Frustum: the near or the far parameter is negative.");
              return newIdentityMatrix();
          }
      //we actually use (near-far) in the matrix below,
      //but for clarity we use far-near, and use minus
          var sx = 2 * near / lx;
          var sy = 2 * near / ly;
          var a = (right + left) / lx;
          var b = (top + bottom) / ly;
          var c1 = -2 * near * far / lz;
          var c2 = -(far + near) / lz;
      ret = [
              sx,  0,  a,  0,
              0 , sy,  b,  0,
              0 ,  0, c2, c1,
              0 ,  0, -1,  0,
      ]
      return ret;
    }
  
    function projectionMatrixOrtho(left, right, bottom, top, near, far, ret) {
      ret = ret || newMatrix();
      var lx = right - left; 
      var ly = top - bottom;
      var lz = far - near;
      ret = [
        2/lx,    0,     0, -(right+left)/lx,
        0,    2/ly,     0, -(top+bottom)/ly,
        0,       0, -2/lz,   -(far+near)/lz,
        0,       0,     0,                1,
      ]
      return ret;
    }
  
      function transpose(src, ret) {
          ret = ret || newMatrix();
          ret = ret || newMatrix();
          var pos0, pos1;
      for(var r = 0; r < matrixSize; ++r)
        for(var c = 0; c < matrixSize; ++c) {
                  pos0 = c * matrixSize + r;
                  pos1 = r * matrixSize + c;
                  ret[pos0] = src[pos1];
              }
          return ret;
    }
      
    function rotationMatrix(angleInRadians, x, y, z, ret) {
      ret = ret || newMatrix();  
      var m = x**2 + y**2 + z**2;
      if(m == 0) return newIdentityMatrix()
      var ux = Math.sqrt(x**2 / m);
      var uy = Math.sqrt(y**2 / m);
      var uz = Math.sqrt(z**2 / m);
      var s = Math.sin(angleInRadians);
      var c = Math.cos(angleInRadians);
      var mc = 1 - c;
      ret = [
                   c + ux**2 * mc, ux * uy * mc -    uz *  s, ux * uz * mc +    uy *  s, 0,
        uy * ux * mc + uz    *  s,            c + uy**2 * mc, uy * uz * mc -    ux *  s, 0,
        uz * ux * mc - uy    *  s, uz * uy * mc +    ux *  s,            c + uz**2 * mc, 0,
        0                        ,                     0,                              0, 1,
      ];
      return ret;
    }
  
      function translationMatrix(x, y, z, ret) {
          ret = ret || newIdentityMatrix();
          var pos;
          var tr =  newVector();
          tr[0] = x;
          tr[1] = y;
          tr[2] = z;
          tr[3] = 1;
          for(var r = 0; r < matrixSize; ++r) {
              pos = r * matrixSize + 3;
              ret[pos] = tr[r];
          }
          return ret;
      }
      
    return {
      matrixSize: matrixSize,
      newMatrix: newMatrix,
      newIdentityMatrix: newIdentityMatrix,
      add: add,
      mulMatrices: mulMatrices,
      mulMatrixVector: mulMatrixVector,
      mulVectorMatrix: mulVectorMatrix,
      mulVectorsInner: mulVectorsInner,
      mul: mul,
      projectionMatrixFrustum: projectionMatrixFrustum,
      projectionMatrixOrtho: projectionMatrixOrtho,
      rotationMatrix: rotationMatrix,
          translationMatrix: translationMatrix,
          transpose: transpose,
    }
  }));
