(this["webpackJsonppse-cime"] = this["webpackJsonppse-cime"] || []).push([
  [0],
  {
    305: function (e, t, n) {},
    306: function (e, t, n) {},
    336: function (e, t, n) {},
    342: function (e, t, n) {
      "use strict";
      n.r(t);
      var r = n(1),
        i = n.n(r),
        a = n(60),
        o = n.n(a),
        c = (n(305), n(11)),
        s = (n(306), n(22)),
        l = n(94),
        u = n(95),
        d = n(202),
        h = n(201),
        f = n(166),
        p = n(14),
        g = (n(336), n(374)),
        j = n(271),
        b = n(259),
        m = n(197),
        O = n(278),
        v = n(265),
        x = n(277),
        S = n(72),
        _ = n(232),
        y = n(84),
        C = n(266),
        k = n(227),
        w = n(230),
        P = n(229),
        L = n(272),
        D = n(117),
        I = n(274),
        R = n(64),
        T = n(164),
        E = n(288),
        M = n.n(E),
        A = n(285),
        N = n.n(A),
        F = n(204),
        V = n.n(F),
        W = n(287),
        z = n.n(W),
        B = n(373),
        q = n(284),
        K = n(286),
        G = n.n(K),
        H = n(200),
        J = n.n(H),
        Y = n(28),
        U = n.n(Y),
        Q = n(59),
        X = "omit",
        Z = "https://cime.caleydoapp.org",
        $ = {},
        ee = {};
      function te(e) {
        var t = arguments.length > 1 && void 0 !== arguments[1] && arguments[1];
        return t ? ee[e] : $[e];
      }
      function ne(e) {
        var t = arguments.length > 1 && void 0 !== arguments[1] && arguments[1],
          n = arguments.length > 2 ? arguments[2] : void 0;
        t ? (ee[e] = n) : ($[e] = n);
      }
      function re(e) {
        return ie.apply(this, arguments);
      }
      function ie() {
        return (ie = Object(Q.a)(
          U.a.mark(function e(t) {
            return U.a.wrap(function (e) {
              for (;;)
                switch ((e.prev = e.next)) {
                  case 0:
                    return e.abrupt("return", t);
                  case 1:
                  case "end":
                    return e.stop();
                }
            }, e);
          })
        )).apply(this, arguments);
      }
      var ae = {};
      function oe(e) {
        return ae[e] ? Object.assign(ae[e]) : null;
      }
      function ce(e, t) {
        ae[e] = t;
      }
      function se(e) {
        if (!e.ok) throw Error(e.statusText);
        return e;
      }
      function le(e) {
        return Object.keys(e).includes("error") && alert(e.error), e;
      }
      function ue() {
        return (ue = Object(Q.a)(
          U.a.mark(function e(t) {
            var n;
            return U.a.wrap(function (e) {
              for (;;)
                switch ((e.prev = e.next)) {
                  case 0:
                    return (
                      (n = Z + "/delete_file/" + t),
                      e.abrupt(
                        "return",
                        fetch(n, { method: "GET", credentials: X })
                          .then(se)
                          .then(function (e) {
                            return e.json();
                          })
                          .then(le)
                          .catch(function (e) {
                            alert(
                              "file could not be deleted. please, try again"
                            ),
                              console.log(e);
                          })
                      )
                    );
                  case 2:
                  case "end":
                    return e.stop();
                }
            }, e);
          })
        )).apply(this, arguments);
      }
      function de() {
        return (de = Object(Q.a)(
          U.a.mark(function e() {
            return U.a.wrap(function (e) {
              for (;;)
                switch ((e.prev = e.next)) {
                  case 0:
                    return (
                      Z + "/get_uploaded_files_list",
                      e.abrupt(
                        "return",
                        fetch(
                          "https://cime.caleydoapp.org/get_uploaded_files_list",
                          { method: "GET", credentials: X }
                        )
                          .then(se)
                          .then(function (e) {
                            return e.json();
                          })
                          .then(le)
                          .catch(function (e) {
                            console.log(e);
                          })
                      )
                    );
                  case 2:
                  case "end":
                    return e.stop();
                }
            }, e);
          })
        )).apply(this, arguments);
      }
      function he() {
        return (he = Object(Q.a)(
          U.a.mark(function e(t, n, r) {
            var i, a, o;
            return U.a.wrap(function (e) {
              for (;;)
                switch ((e.prev = e.next)) {
                  case 0:
                    return (
                      (i = new FormData()).append("smilesA", t),
                      i.append("smilesB", n),
                      (a = Z + "/get_difference_highlight"),
                      (o = r
                        ? fetch(a, {
                            method: "POST",
                            body: i,
                            credentials: X,
                            signal: r.signal,
                          })
                        : fetch(a, {
                            method: "POST",
                            body: i,
                            credentials: X,
                          })),
                      e.abrupt(
                        "return",
                        o
                          .then(se)
                          .then(function (e) {
                            return e.json();
                          })
                          .then(le)
                          .then(function (e) {
                            return console.log(e), e.data;
                          })
                          .catch(function (e) {
                            console.log(e);
                          })
                      )
                    );
                  case 6:
                  case "end":
                    return e.stop();
                }
            }, e);
          })
        )).apply(this, arguments);
      }
      function fe() {
        return (
          (fe = Object(Q.a)(
            U.a.mark(function e(t) {
              var n,
                r,
                i,
                a,
                o,
                c,
                s = arguments;
              return U.a.wrap(function (e) {
                for (;;)
                  switch ((e.prev = e.next)) {
                    case 0:
                      if (
                        ((n = s.length > 1 && void 0 !== s[1] && s[1]),
                        (r = s.length > 2 ? s[2] : void 0),
                        !(i = te(t, n)))
                      ) {
                        e.next = 5;
                        break;
                      }
                      return e.abrupt("return", re(i));
                    case 5:
                      return (
                        (a = new FormData()).append("smiles", t),
                        localStorage.getItem("id") &&
                          a.append("filename", localStorage.getItem("id")),
                        (o = Z + "/get_mol_img"),
                        n && (o += "/highlight"),
                        (c = r
                          ? fetch(o, {
                              method: "POST",
                              body: a,
                              credentials: X,
                              signal: r.signal,
                            })
                          : fetch(o, {
                              method: "POST",
                              body: a,
                              credentials: X,
                            })),
                        e.abrupt(
                          "return",
                          c
                            .then(se)
                            .then(function (e) {
                              return e.json();
                            })
                            .then(le)
                            .then(function (e) {
                              return ne(t, n, e.data), e.data;
                            })
                            .catch(function (e) {
                              console.log(e);
                            })
                        )
                      );
                    case 12:
                    case "end":
                      return e.stop();
                  }
              }, e);
            })
          )),
          fe.apply(this, arguments)
        );
      }
      function pe(e, t) {
        return ge.apply(this, arguments);
      }
      function ge() {
        return (ge = Object(Q.a)(
          U.a.mark(function e(t, n) {
            var r;
            return U.a.wrap(function (e) {
              for (;;)
                switch ((e.prev = e.next)) {
                  case 0:
                    return (
                      localStorage.getItem("id") &&
                        t.append("filename", localStorage.getItem("id")),
                      (r = n
                        ? fetch(Z + "/get_mol_imgs", {
                            method: "POST",
                            body: t,
                            credentials: X,
                            signal:
                              null === n || void 0 === n ? void 0 : n.signal,
                          })
                        : fetch(Z + "/get_mol_imgs", {
                            method: "POST",
                            body: t,
                            credentials: X,
                          })),
                      e.abrupt(
                        "return",
                        r
                          .then(se)
                          .then(function (e) {
                            return e.json();
                          })
                          .then(le)
                          .then(function (e) {
                            return (
                              e.error_smiles.length > 0 &&
                                alert(
                                  "following smiles couldn not be parsed: " +
                                    e.error_smiles
                                ),
                              e
                            );
                          })
                          .catch(function (e) {
                            "AbortError" === e.name
                              ? console.log("Fetch aborted")
                              : (alert("could not load structures"),
                                console.log(e));
                          })
                      )
                    );
                  case 3:
                  case "end":
                    return e.stop();
                }
            }, e);
          })
        )).apply(this, arguments);
      }
      function je() {
        return (je = Object(Q.a)(
          U.a.mark(function e(t, n) {
            var r;
            return U.a.wrap(function (e) {
              for (;;)
                switch ((e.prev = e.next)) {
                  case 0:
                    return (
                      (r = n
                        ? fetch(Z + "/get_common_mol_img", {
                            method: "POST",
                            body: t,
                            signal:
                              null === n || void 0 === n ? void 0 : n.signal,
                          })
                        : fetch(Z + "/get_common_mol_img", {
                            method: "POST",
                            body: t,
                          })),
                      e.abrupt(
                        "return",
                        r
                          .then(se)
                          .then(function (e) {
                            return e.json();
                          })
                          .then(le)
                          .then(function (e) {
                            return e.data;
                          })
                          .catch(function (e) {
                            console.log(e);
                          })
                      )
                    );
                  case 2:
                  case "end":
                    return e.stop();
                }
            }, e);
          })
        )).apply(this, arguments);
      }
      function be(e, t) {
        return me.apply(this, arguments);
      }
      function me() {
        return (me = Object(Q.a)(
          U.a.mark(function e(t, n) {
            var r, i;
            return U.a.wrap(function (e) {
              for (;;)
                switch ((e.prev = e.next)) {
                  case 0:
                    return (
                      (r = new FormData()).append("myFile", t),
                      r.append("file_size", t.size),
                      (i = fetch(Z + "/upload_sdf", {
                        method: "POST",
                        body: r,
                        credentials: X,
                        signal: null === n || void 0 === n ? void 0 : n.signal,
                      })
                        .then(se)
                        .then(function (e) {
                          return e.json();
                        })
                        .then(le)
                        .then(function (e) {
                          localStorage.setItem("id", e.id);
                        })
                        .catch(function (e) {
                          "AbortError" === e.name
                            ? console.log("Fetch aborted")
                            : (alert(
                                "error when uploading file. it might be too big"
                              ),
                              console.log(e));
                        })),
                      e.abrupt("return", i)
                    );
                  case 5:
                  case "end":
                    return e.stop();
                }
            }, e);
          })
        )).apply(this, arguments);
      }
      function Oe() {
        return ve.apply(this, arguments);
      }
      function ve() {
        return (
          (ve = Object(Q.a)(
            U.a.mark(function e() {
              var t,
                n,
                r,
                i,
                a,
                o,
                c = arguments;
              return U.a.wrap(function (e) {
                for (;;)
                  switch ((e.prev = e.next)) {
                    case 0:
                      if (
                        ((t = c.length > 0 && void 0 !== c[0] && c[0]),
                        (n = c.length > 1 && void 0 !== c[1] ? c[1] : ""),
                        (r = c.length > 2 ? c[2] : void 0),
                        t)
                      ) {
                        e.next = 7;
                        break;
                      }
                      if (
                        !(
                          (i = oe("representation_list_" + n)) &&
                          i.rep_list.length > 0
                        )
                      ) {
                        e.next = 7;
                        break;
                      }
                      return e.abrupt("return", re(i));
                    case 7:
                      return (
                        (a = Z + "/get_atom_rep_list"),
                        localStorage.getItem("id") &&
                          (a += "/" + localStorage.getItem("id")),
                        (o = r
                          ? fetch(a, {
                              method: "GET",
                              credentials: X,
                              signal: r.signal,
                            })
                          : fetch(a, { method: "GET", credentials: X })),
                        e.abrupt(
                          "return",
                          o
                            .then(se)
                            .then(function (e) {
                              return e.json();
                            })
                            .then(le)
                            .then(function (e) {
                              return ce("representation_list_" + n, e), e;
                            })
                            .catch(function (e) {
                              console.log(e);
                            })
                        )
                      );
                    case 11:
                    case "end":
                      return e.stop();
                  }
              }, e);
            })
          )),
          ve.apply(this, arguments)
        );
      }
      var xe = n(2),
        Se = Object(T.connect)(
          function (e) {
            var t;
            return {
              dataset: e.dataset,
              hoverSettings: e.hoverSettings,
              rdkitSettings: e.rdkitSettings,
              columns:
                null === (t = e.dataset) || void 0 === t ? void 0 : t.columns,
            };
          },
          function (e) {
            return {
              setCurrentAggregation: function (t) {
                return e(Object(s.selectVectors)(t, !1));
              },
              setHoverstate: function (t, n) {
                return e(Object(s.setHoverState)(t, n));
              },
            };
          }
        ),
        _e = Se(function (e) {
          var t = Object(s.useCancellablePromise)(),
            n = t.cancellablePromise,
            i = t.cancelPromises;
          if (e.mcs_only) {
            var a = r.useState(
                Object(xe.jsx)("div", { children: "loading..." })
              ),
              o = Object(c.a)(a, 2),
              l = o[0],
              u = o[1],
              d = Pe(e.columns);
            return (
              r.useEffect(
                function () {
                  if ((i(), d in e.columns)) {
                    var t,
                      r = new AbortController();
                    if (e.diff && e.selection_ref) {
                      t = (function (e, t, n) {
                        return he.apply(this, arguments);
                      })(
                        e.selection.map(function (e) {
                          return e[d];
                        }),
                        e.selection_ref.map(function (e) {
                          return e[d];
                        }),
                        r
                      );
                    } else {
                      var a = new FormData();
                      e.selection.every(function (e) {
                        return a.append("smiles_list", e[d]), !0;
                      }),
                        (t = (function (e, t) {
                          return je.apply(this, arguments);
                        })(a, r));
                    }
                    Object(R.trackPromise)(
                      n(
                        t.then(function (e) {
                          e.length > 100
                            ? u(function () {
                                return Object(xe.jsx)("div", {
                                  style: {
                                    width: 200,
                                    height: 200,
                                    backgroundSize: "contain",
                                    backgroundPosition: "center",
                                    backgroundRepeat: "no-repeat",
                                    backgroundImage:
                                      "url('data:image/jpg;base64,".concat(
                                        e,
                                        "')"
                                      ),
                                  },
                                });
                              })
                            : u(function () {
                                return Object(xe.jsx)("div", { children: e });
                              });
                        }),
                        r
                      )
                    );
                  }
                },
                [e.selection, e.selection_ref, e.mcs_only]
              ),
              d in e.columns
                ? Object(xe.jsx)("div", { children: l })
                : Object(xe.jsx)("div", { children: "No SMILES column found." })
            );
          }
          var h = r.useState(!1),
            f = Object(c.a)(h, 2),
            m = f[0],
            O = f[1],
            v = r.useState(["Common Substructure"]),
            x = Object(c.a)(v, 2),
            S = x[0],
            _ = x[1],
            y = r.useState([0]),
            C = Object(c.a)(y, 2),
            k = C[0],
            w = C[1],
            P = r.useState(["Common Substructure"]),
            L = Object(c.a)(P, 2),
            D = L[0],
            I = L[1],
            T = function () {
              var t =
                arguments.length > 0 && void 0 !== arguments[0] && arguments[0];
              if (t || S.length <= 1) {
                var r = "global_loading_indicator",
                  i = new AbortController();
                Object(R.trackPromise)(
                  n(
                    Oe(t, e.dataset.info.path, i).then(function (e) {
                      if (e.rep_list.length > 0) {
                        var t = Object(p.a)(e.rep_list);
                        t.splice(0, 0, "Common Substructure"), _(t);
                      }
                    }),
                    i
                  ),
                  r
                );
              }
            };
          r.useEffect(function () {
            i(), e.aggregate && T();
          }, []);
          var E = function (e, t) {
              if (S.includes(e)) {
                var n = Object(p.a)(D);
                (n[k.indexOf(t)] = e), I(n);
              }
            },
            M = r.useRef(),
            A = r.useRef();
          return e.aggregate
            ? Object(xe.jsxs)(g.a, {
                className: "ParentChem",
                paddingBottom: 3,
                children: [
                  e.aggregate &&
                    Object(xe.jsxs)(g.a, {
                      paddingLeft: 2,
                      paddingRight: 2,
                      children: [
                        Object(xe.jsx)(j.a, {
                          title: "Summary Settings",
                          children: Object(xe.jsxs)(b.a, {
                            style: { color: "gray" },
                            ref: M,
                            onClick: function () {
                              return O(!0);
                            },
                            children: [
                              Object(xe.jsx)(N.a, {}),
                              "\xa0 Settings",
                            ],
                          }),
                        }),
                        Object(xe.jsx)(Ie, {
                          open: m,
                          setOpen: O,
                          anchorEl: M.current,
                          refreshRepList: function () {
                            T(!0);
                          },
                        }),
                        Object(xe.jsx)(j.a, {
                          title: "Add Component",
                          children: Object(xe.jsxs)(b.a, {
                            style: { color: "gray" },
                            onClick: function () {
                              return (function () {
                                var e = Object(p.a)(k);
                                e.push(
                                  Math.max.apply(Math, Object(p.a)(e)) + 1
                                ),
                                  w(e);
                                var t = Object(p.a)(D);
                                t.push("Common Substructure"), I(t);
                              })();
                            },
                            children: [
                              Object(xe.jsx)(G.a, {}),
                              "\xa0 Add View",
                            ],
                          }),
                        }),
                      ],
                    }),
                  Object(xe.jsxs)("div", {
                    ref: A,
                    className: "chemComponents",
                    children: [
                      k.length > 1 &&
                        Object(xe.jsx)("div", {
                          style: {
                            width: (e.rdkitSettings.width + 20) * k.length,
                          },
                          children: k.map(function (t, n) {
                            return Object(xe.jsx)(
                              "div",
                              {
                                style: {
                                  width: e.rdkitSettings.width + 20,
                                  float: "left",
                                },
                                children: Object(xe.jsx)(ke, {
                                  chemRef: A,
                                  setCurrentRep: function (e) {
                                    return E(e, t);
                                  },
                                  currentRep: D[n],
                                  removeComponent: function () {
                                    return (function (e) {
                                      var t = Object(p.a)(k),
                                        n = Object(p.a)(D),
                                        r = t.indexOf(e);
                                      r > -1 &&
                                        (t.splice(r, 1), n.splice(r, 1)),
                                        w(t),
                                        I(n);
                                    })(t);
                                  },
                                  id: t,
                                  rep_list: S,
                                  selection: e.selection,
                                  aggregate: e.aggregate,
                                }),
                              },
                              t
                            );
                          }),
                        }),
                      k.length <= 1 &&
                        Object(xe.jsx)("div", {
                          children: Object(xe.jsx)(
                            "div",
                            {
                              style: { minWidth: e.rdkitSettings.width },
                              children: Object(xe.jsx)(ke, {
                                chemRef: A,
                                setCurrentRep: function (e) {
                                  return E(e, k[0]);
                                },
                                currentRep: D[0],
                                id: k[0],
                                rep_list: S,
                                selection: e.selection,
                                aggregate: e.aggregate,
                              }),
                            },
                            k[0]
                          ),
                        }),
                    ],
                  }),
                ],
              })
            : Object(xe.jsx)(ke, {
                id: -1,
                rep_list: S,
                selection: e.selection,
                aggregate: e.aggregate,
              });
        }),
        ye = "chemlegend_loading_area",
        Ce = "chemdetail",
        ke = Se(
          (function (e) {
            Object(d.a)(n, e);
            var t = Object(h.a)(n);
            function n(e) {
              var r;
              return (
                Object(l.a)(this, n),
                ((r = t.call(this, e)).anchorRef = void 0),
                (r.state = { checkedList: [] }),
                r
              );
            }
            return (
              Object(u.a)(n, [
                {
                  key: "render",
                  value: function () {
                    var e = this,
                      t = function (t) {
                        var n = Object(q.isFunction)(t)
                          ? t(e.state.checkedList)
                          : t;
                        e.setState(
                          Object(f.a)(
                            Object(f.a)({}, e.state),
                            {},
                            { checkedList: n }
                          )
                        );
                      };
                    return this.props.aggregate
                      ? Object(xe.jsxs)("div", {
                          className: "ParentImg",
                          children: [
                            Object(xe.jsx)(g.a, {
                              paddingLeft: 2,
                              paddingTop: 1,
                              paddingRight: 2,
                              children: Object(xe.jsx)(Re, {
                                value: this.props.currentRep,
                                onChange: this.props.setCurrentRep,
                                rep_list: this.props.rep_list,
                                hoverSettings: this.props.hoverSettings,
                              }),
                            }),
                            Object(xe.jsxs)(g.a, {
                              paddingLeft: 2,
                              paddingTop: 1,
                              paddingRight: 2,
                              children: [
                                Object(xe.jsxs)(b.a, {
                                  size: "small",
                                  variant: "outlined",
                                  onClick: function () {
                                    !(function () {
                                      var n = e.props.selection.filter(
                                        function (t, n) {
                                          return e.state.checkedList[n];
                                        }
                                      );
                                      n.length > 0
                                        ? (t([]),
                                          e.props.setCurrentAggregation(
                                            n.map(function (e) {
                                              return e.__meta__.meshIndex;
                                            })
                                          ))
                                        : alert(
                                            "Please, select at least one Compound in the Summary View to filter."
                                          );
                                    })();
                                  },
                                  children: [
                                    Object(xe.jsx)(z.a, { fontSize: "small" }),
                                    "\xa0Confirm Selection",
                                  ],
                                }),
                                this.props.removeComponent &&
                                  Object(xe.jsx)(m.a, {
                                    onClick: this.props.removeComponent,
                                    children: Object(xe.jsx)(J.a, {}),
                                  }),
                              ],
                            }),
                            Object(xe.jsx)(s.LoadingIndicatorView, {
                              area: ye + this.props.id,
                            }),
                            Object(xe.jsx)(De, {
                              chemRef: this.props.chemRef,
                              id: this.props.id,
                              setCheckedList: t,
                              selection: this.props.selection,
                              columns: this.props.columns,
                              aggregate: this.props.aggregate,
                              current_rep: this.props.currentRep,
                              handleMouseEnter: function (t) {
                                var n = null;
                                t >= 0 && (n = e.props.selection[t]),
                                  e.props.setHoverstate(n, Ce);
                              },
                              handleMouseOut: function () {
                                e.props.setHoverstate(null, Ce);
                              },
                            }),
                          ],
                        })
                      : Object(xe.jsx)("div", {
                          children: Object(xe.jsx)(De, {
                            id: this.props.id,
                            selection: this.props.selection,
                            columns: this.props.columns,
                            aggregate: this.props.aggregate,
                          }),
                        });
                  },
                },
              ]),
              n
            );
          })(r.Component)
        );
      function we(e, t, n, r, i, a) {
        var o = Pe(e.columns);
        if (o in e.columns)
          if ((t(Object(xe.jsx)("div", {})), e.selection.length > 0))
            if (e.aggregate) {
              var c = new FormData();
              c.append("current_rep", e.current_rep),
                e.selection.forEach(function (e) {
                  c.append("smiles_list", e[o]);
                }),
                c.append("contourLines", e.rdkitSettings.contourLines),
                c.append("scale", e.rdkitSettings.scale),
                c.append("sigma", e.rdkitSettings.sigma),
                c.append("showMCS", e.rdkitSettings.showMCS),
                c.append("width", e.rdkitSettings.width),
                c.append("doAlignment", e.rdkitSettings.doAlignment);
              var s = new AbortController();
              Object(R.trackPromise)(
                i(pe(c, s), s).then(function (i) {
                  var o = i.img_lst.map(function (t, i) {
                    return (
                      a(function (e) {
                        var t = Object(p.a)(e);
                        return t.length <= i && t.push(!1), t;
                      }),
                      Object(xe.jsxs)(
                        O.a,
                        {
                          className: "legend_multiple",
                          item: !0,
                          children: [
                            Object(xe.jsx)(v.a, {
                              labelPlacement: "bottom",
                              control: Object(xe.jsx)(x.a, {
                                color: "primary",
                                onChange: function (e) {
                                  !(function (e, t) {
                                    a(function (n) {
                                      return n.map(function (n, r) {
                                        return r === e ? t : n;
                                      });
                                    });
                                  })(i, e.target.checked);
                                },
                              }),
                              label: Object(xe.jsx)("img", {
                                src: "data:image/jpeg;base64," + t,
                                onMouseEnter: function () {
                                  n(i);
                                },
                                onMouseOver: function () {
                                  n(i);
                                },
                                onMouseLeave: function () {
                                  r();
                                },
                              }),
                            }),
                            Object(xe.jsxs)(S.a, {
                              style: { paddingLeft: 5 },
                              variant: "subtitle2",
                              children: ["ID: ", e.selection[i].ID],
                            }),
                          ],
                        },
                        i
                      )
                    );
                  });
                  t(o);
                }),
                ye + e.id
              );
            } else {
              var l = e.selection[0],
                u = new AbortController();
              i(
                (function (e) {
                  return fe.apply(this, arguments);
                })(l[o], !1, u),
                u
              )
                .then(function (e) {
                  t(
                    Object(xe.jsxs)("div", {
                      children: [
                        Object(xe.jsx)("img", {
                          className: "legend_single",
                          src: "data:image/jpeg;base64," + e,
                        }),
                        Object(xe.jsxs)(S.a, {
                          style: { paddingLeft: 5 },
                          variant: "subtitle2",
                          children: ["ID: ", l.ID],
                        }),
                      ],
                    })
                  );
                })
                .catch(function (e) {
                  return console.log(e);
                });
            }
          else t(Object(xe.jsx)("div", { children: "No Selection" }));
        else t(Object(xe.jsx)("div", { children: "No SMILES column found" }));
      }
      function Pe(e) {
        var t = "SMILES";
        if (!(t in e)) {
          var n = Object.keys(e);
          for (var r in n) {
            var i = n[r];
            if (i.toLowerCase().includes("smiles")) {
              t = i;
              break;
            }
          }
        }
        return t;
      }
      function Le(e) {
        e && e.style && (e.style.border = "solid white 1px");
      }
      var De = Object(T.connect)(
        function (e) {
          return { hoverState: e.hoverState, rdkitSettings: e.rdkitSettings };
        },
        function (e) {
          return {};
        }
      )(function (e) {
        var t = e.chemRef,
          n = e.id,
          i = e.hoverState,
          a = e.selection,
          o = e.columns,
          l = e.aggregate,
          u = e.handleMouseEnter,
          d = e.handleMouseOut,
          h = e.current_rep,
          f = e.setCheckedList,
          p = e.rdkitSettings,
          g = r.useState(Object(xe.jsx)("div", {})),
          j = Object(c.a)(g, 2),
          b = j[0],
          m = j[1],
          v = r.useRef(),
          x = Object(s.useCancellablePromise)(),
          S = x.cancellablePromise,
          _ = x.cancelPromises;
        return (
          r.useEffect(
            function () {
              _();
            },
            [a, h]
          ),
          r.useEffect(
            function () {
              f && f([]),
                we(
                  {
                    id: n,
                    columns: o,
                    aggregate: l,
                    current_rep: h,
                    selection: a,
                    rdkitSettings: p,
                  },
                  m,
                  u,
                  d,
                  S,
                  f
                );
            },
            [a]
          ),
          r.useEffect(
            function () {
              l &&
                (function (e, t) {
                  var n = Pe(e.columns);
                  if (n in e.columns) {
                    var r = e.imgContainer.childNodes;
                    if (e.selection.length == r.length) {
                      e.imgContainer.style.display = "none";
                      var i = new FormData();
                      i.append("current_rep", e.current_rep),
                        e.selection.forEach(function (e) {
                          i.append("smiles_list", e[n]);
                        }),
                        i.append("contourLines", e.rdkitSettings.contourLines),
                        i.append("scale", e.rdkitSettings.scale),
                        i.append("sigma", e.rdkitSettings.sigma),
                        i.append("showMCS", e.rdkitSettings.showMCS),
                        i.append("width", e.rdkitSettings.width),
                        i.append("doAlignment", e.rdkitSettings.doAlignment);
                      var a = new AbortController();
                      Object(R.trackPromise)(
                        t(pe(i, a), a).then(function (t) {
                          t.img_lst.map(function (e, t) {
                            r[t].getElementsByTagName("img")[0].src =
                              "data:image/jpeg;base64," + e;
                          }),
                            (e.imgContainer.style.display = "flex");
                        }),
                        ye + e.id
                      );
                    }
                  }
                })(
                  {
                    id: n,
                    columns: o,
                    current_rep: h,
                    selection: a,
                    imgContainer:
                      null === v || void 0 === v ? void 0 : v.current,
                    rdkitSettings: p,
                  },
                  S
                );
            },
            [h, p.refresh]
          ),
          r.useEffect(
            function () {
              if (l) {
                var e = null === t || void 0 === t ? void 0 : t.current,
                  n = e.getElementsByClassName("chem-grid")[0].childNodes;
                if (i && i.data) {
                  var r = a.findIndex(function (e) {
                    return (
                      e &&
                      e.__meta__ &&
                      i.data.__meta__ &&
                      e.__meta__.meshIndex == i.data.__meta__.meshIndex
                    );
                  });
                  if (r >= 0 && n.length > 0) {
                    for (var o in n) {
                      Le(n[o]);
                    }
                    (s = n[r]) &&
                      s.style &&
                      (s.style.border = "solid black 1px"),
                      i.updater != Ce &&
                        e &&
                        n[r] &&
                        (e.scrollTop = n[r].offsetTop - e.offsetTop);
                  }
                } else
                  for (var c in n) {
                    Le(n[c]);
                  }
              }
              var s;
            },
            [i.data, i.updater]
          ),
          Object(xe.jsx)("div", {
            className: "chemContainer",
            children: Object(xe.jsx)(O.a, {
              ref: v,
              className: "chem-grid",
              container: !0,
              children: b,
            }),
          })
        );
      });
      var Ie = Object(T.connect)(
          function (e) {
            return { rdkitSettings: e.rdkitSettings };
          },
          function (e) {
            return {
              setContourLines: function (t) {
                return e(Object(s.setRDKit_contourLines)(t));
              },
              setScale: function (t) {
                return e(Object(s.setRDKit_scale)(t));
              },
              setSigma: function (t) {
                return e(Object(s.setRDKit_sigma)(t));
              },
              setShowMCS: function (t) {
                return e(Object(s.setRDKit_showMCS)(t));
              },
              setWidth: function (t) {
                return e(Object(s.setRDKit_width)(t));
              },
              setRefresh: function (t) {
                return e(Object(s.setRDKit_refresh)(t));
              },
              setDoAlignment: function (t) {
                return e(Object(s.setRDKit_doAlignment)(t));
              },
            };
          }
        )(function (e) {
          var t = e.open,
            n = e.setOpen,
            r = e.anchorEl,
            i = e.refreshRepList,
            a = e.rdkitSettings,
            o = e.setContourLines,
            c = e.setScale,
            s = e.setSigma,
            l = e.setShowMCS,
            u = (e.setWidth, e.setRefresh),
            d = e.setDoAlignment;
          return Object(xe.jsx)(_.a, {
            disablePortal: !0,
            id: "dialog to open",
            open: t,
            anchorEl: r,
            onClose: function () {
              return n(function () {
                return !1;
              });
            },
            anchorOrigin: { vertical: "bottom", horizontal: "right" },
            transformOrigin: { vertical: "bottom", horizontal: "left" },
            children: Object(xe.jsx)("div", {
              children: Object(xe.jsx)(y.a, {
                style: { padding: 10, minWidth: 300 },
                children: Object(xe.jsxs)(C.a, {
                  children: [
                    Object(xe.jsxs)(b.a, {
                      size: "small",
                      variant: "outlined",
                      "aria-label": "Refresh Representation List",
                      onClick: function () {
                        return i(!0);
                      },
                      children: [
                        Object(xe.jsx)(M.a, {}),
                        "Refresh Representation List",
                      ],
                    }),
                    Object(xe.jsx)(S.a, {
                      variant: "subtitle2",
                      gutterBottom: !0,
                      children: "RDKit Settings",
                    }),
                    Object(xe.jsxs)(k.a, {
                      children: [
                        Object(xe.jsxs)(w.a, {
                          shrink: !0,
                          htmlFor: "contourLinesInput",
                          children: [
                            "Contour Lines ",
                            Object(xe.jsx)(j.a, {
                              title:
                                "Number of Contour Lines [0; \u221e] \u2208 \u2115",
                              children: Object(xe.jsx)(V.a, {
                                fontSize: "small",
                              }),
                            }),
                          ],
                        }),
                        Object(xe.jsx)(P.a, {
                          id: "contourLinesInput",
                          type: "number",
                          value: a.contourLines,
                          onChange: function (e) {
                            var t = parseInt(e.target.value);
                            isNaN(t) ? o(e.target.value) : o(Math.max(t, 0));
                          },
                        }),
                      ],
                    }),
                    Object(xe.jsxs)(k.a, {
                      children: [
                        Object(xe.jsxs)(w.a, {
                          shrink: !0,
                          htmlFor: "ScaleInput",
                          children: [
                            "Scale ",
                            Object(xe.jsx)(j.a, {
                              title: "Weight Scale [-1; \u221e] \u2208 \u211d",
                              children: Object(xe.jsx)(V.a, {
                                fontSize: "small",
                              }),
                            }),
                          ],
                        }),
                        Object(xe.jsx)(P.a, {
                          id: "ScaleInput",
                          type: "number",
                          value: a.scale,
                          onChange: function (e) {
                            var t = parseFloat(e.target.value);
                            isNaN(t) ? c(e.target.value) : c(Math.max(t, -1));
                          },
                        }),
                      ],
                    }),
                    Object(xe.jsxs)(k.a, {
                      children: [
                        Object(xe.jsxs)(w.a, {
                          shrink: !0,
                          htmlFor: "SigmaInput",
                          children: [
                            "Sigma ",
                            Object(xe.jsx)(j.a, {
                              title:
                                "Sigma for Gaussian ]0; \u221e] \u2208 \u211d. Default of 0 signals the algorithm to infer the value.",
                              children: Object(xe.jsx)(V.a, {
                                fontSize: "small",
                              }),
                            }),
                          ],
                        }),
                        Object(xe.jsx)(P.a, {
                          id: "SigmaInput",
                          type: "number",
                          inputProps: { step: 0.1 },
                          value: a.sigma,
                          onChange: function (e) {
                            var t = parseFloat(e.target.value);
                            isNaN(t) ? s(e.target.value) : s(Math.max(t, 0));
                          },
                        }),
                      ],
                    }),
                    Object(xe.jsx)(v.a, {
                      control: Object(xe.jsx)(L.a, {
                        color: "primary",
                        checked: a.showMCS,
                        onChange: function (e, t) {
                          l(t);
                        },
                      }),
                      label: "Show MCS",
                    }),
                    Object(xe.jsx)(v.a, {
                      control: Object(xe.jsx)(L.a, {
                        color: "primary",
                        checked: a.doAlignment,
                        onChange: function (e, t) {
                          d(t);
                        },
                      }),
                      label: "Align Structure",
                    }),
                    Object(xe.jsx)(S.a, {
                      style: { paddingTop: 10 },
                      gutterBottom: !0,
                      children: "Image Width",
                    }),
                    Object(xe.jsx)(b.a, {
                      style: { marginTop: 3, maxWidth: 150 },
                      size: "small",
                      variant: "outlined",
                      onClick: function () {
                        u((a.refresh += 1));
                      },
                      children: "Apply Settings",
                    }),
                  ],
                }),
              }),
            }),
          });
        }),
        Re = function (e) {
          var t = e.rep_list.map(function (e) {
              var t = e.split("_"),
                n = t.pop(),
                r = t.join("_");
              return {
                group: (r = (r = r.replace("atom.dprop.", "")).replace(
                  "atom.dprop",
                  ""
                )),
                value: e,
                inputValue: n,
              };
            }),
            n = Object(D.a)({
              stringify: function (e) {
                return e.value;
              },
            });
          return Object(xe.jsx)(B.a, {
            size: "small",
            className: e.className,
            filterOptions: n,
            onChange: function (t, n) {
              n && e.onChange(n.value);
            },
            disablePortal: e.hoverSettings.windowMode == s.WindowMode.Extern,
            options: t.sort(function (e, t) {
              return -t.group.localeCompare(e.group);
            }),
            groupBy: function (e) {
              return e.group;
            },
            getOptionLabel: function (e) {
              return e.inputValue;
            },
            style: { maxWidth: 300 },
            defaultValue: t[0],
            renderInput: function (e) {
              return Object(xe.jsx)(
                I.a,
                Object(f.a)(
                  Object(f.a)({}, e),
                  {},
                  { label: "Choose Representation", variant: "outlined" }
                )
              );
            },
          });
        },
        Te = (function (e) {
          Object(d.a)(n, e);
          var t = Object(h.a)(n);
          function n() {
            var e;
            Object(l.a)(this, n);
            for (var r = arguments.length, i = new Array(r), a = 0; a < r; a++)
              i[a] = arguments[a];
            return (
              ((e = t.call.apply(t, [this].concat(i))).type =
                s.DatasetType.Chem),
              e
            );
          }
          return (
            Object(u.a)(n, [
              {
                key: "createFingerprint",
                value: function (e, t, n) {
                  return Object(xe.jsx)(_e, { selection: e, aggregate: n });
                },
              },
            ]),
            n
          );
        })(s.PSEPlugin),
        Ee = n(341);
      var Me = (function () {
          function e() {
            Object(l.a)(this, e),
              (this.vectors = []),
              (this.datasetType = s.DatasetType.None),
              (this.loading_area = "global_loading_indicator");
          }
          return (
            Object(u.a)(e, [
              {
                key: "resolvePath",
                value: function (e, t, n, r, i) {
                  var a = this;
                  e.uploaded
                    ? (localStorage.setItem("id", e.path),
                      this.loadCSV(t, e, n, r, i))
                    : Object(R.trackPromise)(
                        fetch(e.path, {
                          signal:
                            null === i || void 0 === i ? void 0 : i.signal,
                        })
                          .then(function (e) {
                            return e.blob();
                          })
                          .then(function (e) {
                            return a.resolveContent(e, t, n, r, i);
                          })
                          .catch(function (e) {
                            console.log(e);
                          }),
                        this.loading_area
                      );
                },
              },
              {
                key: "resolveContent",
                value: function (e, t, n, r, i) {
                  var a = this,
                    o = n ? n(be(e, i), i) : be(e, i);
                  Object(R.trackPromise)(
                    o
                      .then(function () {
                        a.loadCSV(
                          t,
                          { display: "", type: a.datasetType, path: e.name },
                          n,
                          r,
                          i
                        );
                      })
                      .catch(function (e) {
                        console.log(e);
                      }),
                    this.loading_area
                  );
                },
              },
              {
                key: "loadCSV",
                value: function (e, t, n, r, i) {
                  var a = this,
                    o = Z + "/get_csv/",
                    c = localStorage.getItem("id");
                  void 0 !== c && (o += c), (o += "/"), (o += r);
                  var l = n
                    ? n(
                        Ee.csv(o, {
                          credentials: X,
                          signal:
                            null === i || void 0 === i ? void 0 : i.signal,
                        }),
                        i
                      )
                    : Ee.csv(o, {
                        credentials: X,
                        signal: null === i || void 0 === i ? void 0 : i.signal,
                      });
                  Object(R.trackPromise)(
                    l
                      .then(function (n) {
                        (a.vectors = (function (e) {
                          return e.map(function (e) {
                            return s.AVector.create(e);
                          });
                        })(n)),
                          (a.datasetType = s.DatasetType.Chem),
                          new s.CSVLoader().resolve(
                            e,
                            a.vectors,
                            a.datasetType,
                            t
                          );
                      })
                      .catch(function (e) {
                        console.log(e);
                      }),
                    this.loading_area
                  );
                },
              },
            ]),
            e
          );
        })(),
        Ae = n(260),
        Ne = n(264),
        Fe = n(262),
        Ve = n(263),
        We = n(261);
      function ze(e) {
        var t = e.openSDFDialog,
          n = e.handleClose,
          r = i.a.useState(""),
          a = Object(c.a)(r, 2),
          o = a[0],
          s = a[1];
        return Object(xe.jsxs)(Ae.a, {
          maxWidth: "lg",
          open: t,
          onClose: function () {
            return n(null);
          },
          children: [
            Object(xe.jsx)(Ne.a, { children: "Specify Modifiers" }),
            Object(xe.jsxs)(Fe.a, {
              children: [
                Object(xe.jsxs)(Ve.a, {
                  children: [
                    'Manually specify modifiers separated by semicolons e.g. "pred;fp;latent". ',
                    Object(xe.jsx)("br", {}),
                    "You can also leave this field empty, if the modifiers are included by default. ",
                    Object(xe.jsx)("br", {}),
                    'The following modifiers are included by default: "pred", "predicted", "measured", "fingerprint", "rep".',
                  ],
                }),
                Object(xe.jsx)(I.a, {
                  autoFocus: !0,
                  margin: "dense",
                  id: "modifiers",
                  label: "Modifiers",
                  value: o,
                  onChange: function (e) {
                    s(e.target.value);
                  },
                  fullWidth: !0,
                }),
              ],
            }),
            Object(xe.jsxs)(We.a, {
              children: [
                Object(xe.jsx)(b.a, {
                  onClick: function () {
                    return n(null);
                  },
                  children: "Cancel",
                }),
                Object(xe.jsx)(b.a, {
                  onClick: function () {
                    return n(o);
                  },
                  children: "Start",
                }),
              ],
            }),
          ],
        });
      }
      var Be = function (e) {
          var t = e.onChange,
            n = e.cancellablePromise,
            r = e.abort_controller,
            a = i.a.useState(null),
            o = Object(c.a)(a, 2),
            l = o[0],
            u = o[1],
            d = i.a.useState(!1),
            h = Object(c.a)(d, 2),
            f = h[0],
            p = h[1];
          return Object(xe.jsxs)(O.a, {
            container: !0,
            item: !0,
            alignItems: "stretch",
            justifyContent: "center",
            direction: "column",
            style: { padding: "16px" },
            children: [
              Object(xe.jsx)(s.DragAndDrop, {
                accept: "image/*",
                handleDrop: function (e) {
                  if (!(null == e || e.length <= 0)) {
                    var n = e[0],
                      r = n.name;
                    if (r.endsWith("sdf")) u(n), p(!0);
                    else {
                      var i = new FileReader();
                      (i.onload = function (e) {
                        var n,
                          i =
                            null === e ||
                            void 0 === e ||
                            null === (n = e.target) ||
                            void 0 === n
                              ? void 0
                              : n.result;
                        r.endsWith("json")
                          ? new s.JSONLoader().resolveContent(i, t)
                          : new s.CSVLoader().resolveContent(i, t);
                      }),
                        i.readAsText(n);
                    }
                  }
                },
                children: Object(xe.jsx)("div", { style: { height: 200 } }),
              }),
              Object(xe.jsx)(ze, {
                openSDFDialog: f,
                handleClose: function (e) {
                  p(!1),
                    null !== e &&
                      ((r = new AbortController()),
                      new Me().resolveContent(l, t, n, e, r));
                },
              }),
            ],
          });
        },
        qe = n(267),
        Ke = n(225),
        Ge = n(268),
        He = n(270),
        Je = n(231),
        Ye = n(289),
        Ue = function (e) {
          var t = e.onChange,
            n = e.refresh,
            r = i.a.useState(null),
            a = Object(c.a)(r, 2),
            o = a[0],
            l = a[1],
            u = Object(s.useCancellablePromise)().cancellablePromise;
          i.a.useEffect(
            function () {
              h();
            },
            [n]
          );
          var d = "update_uploaded_files_list";
          function h() {
            Object(R.trackPromise)(
              u(
                (function () {
                  return de.apply(this, arguments);
                })()
              )
                .then(function (e) {
                  e && Object.keys(e).includes("file_list") && l(e.file_list);
                })
                .catch(function (e) {
                  return console.log(e);
                }),
              d
            );
          }
          var f = function (e) {
            u(
              (function (e) {
                return ue.apply(this, arguments);
              })(e)
            )
              .then(function (t) {
                if (t && "true" == t.deleted) {
                  var n = Object(p.a)(o),
                    r = n.indexOf(e);
                  r > -1 && n.splice(r, 1), l(n);
                }
              })
              .catch(function (e) {
                return console.log(e);
              });
          };
          return (
            o &&
            Object(xe.jsx)("div", {
              children: Object(xe.jsxs)(O.a, {
                item: !0,
                style: {
                  overflowY: "auto",
                  flex: "1 1 auto",
                  maxHeight: "400px",
                },
                children: [
                  Object(xe.jsxs)(qe.a, {
                    subheader: Object(xe.jsx)("li", {}),
                    style: { backgroundColor: "white" },
                    children: [
                      !tt &&
                        Object(xe.jsxs)(Ke.a, {
                          children: [
                            "Uploaded Files ",
                            Object(xe.jsx)(b.a, {
                              onClick: function () {
                                return h();
                              },
                              children: Object(xe.jsx)(Ye.a, {}),
                            }),
                          ],
                        }),
                      tt &&
                        Object(xe.jsx)(Ke.a, { children: "Select Dataset" }),
                      o.map(function (e) {
                        return Object(xe.jsxs)(
                          Ge.a,
                          {
                            style: { maxWidth: "270px" },
                            button: !0,
                            onClick: function () {
                              var n;
                              (n = {
                                display: e,
                                path: e,
                                type: s.DatasetType.Chem,
                                uploaded: !0,
                              }),
                                t(n);
                            },
                            children: [
                              Object(xe.jsx)(He.a, { primary: e }),
                              !tt &&
                                Object(xe.jsx)(Je.a, {
                                  onClick: function () {
                                    f(e);
                                  },
                                  children: Object(xe.jsx)(m.a, {
                                    edge: "end",
                                    "aria-label": "delete",
                                    children: Object(xe.jsx)(J.a, {}),
                                  }),
                                }),
                            ],
                          },
                          e
                        );
                      }),
                    ],
                  }),
                  Object(xe.jsx)(s.LoadingIndicatorView, { area: d }),
                ],
              }),
            })
          );
        };
      function Qe(e) {
        var t = e.onDataSelected,
          n = i.a.useState(null),
          r = Object(c.a)(n, 2),
          a = r[0],
          o = r[1],
          l = i.a.useState(!1),
          u = Object(c.a)(l, 2),
          d = u[0],
          h = u[1],
          f = i.a.useState(0),
          p = Object(c.a)(f, 2),
          j = p[0],
          b = p[1],
          m = Object(s.useCancellablePromise)(),
          O = m.cancellablePromise,
          v = m.cancelPromises,
          x = new AbortController();
        return Object(xe.jsxs)("div", {
          style: { display: "flex", flexDirection: "column", height: "100%" },
          children: [
            Object(xe.jsx)(g.a, {
              paddingLeft: 2,
              paddingTop: 2,
              children: Object(xe.jsx)(S.a, {
                variant: "subtitle2",
                gutterBottom: !0,
                children: "Custom Datasets (Drag and Drop)",
              }),
            }),
            Object(xe.jsx)(Be, {
              onChange: function (e, n) {
                t(e, n), b(j + 1);
              },
              cancellablePromise: O,
              abort_controller: x,
            }),
            Object(xe.jsx)(g.a, {
              paddingLeft: 2,
              paddingTop: 2,
              children: Object(xe.jsx)(S.a, {
                variant: "subtitle2",
                gutterBottom: !0,
                children: "Predefined Datasets",
              }),
            }),
            Object(xe.jsx)(Ue, {
              onChange: function (e) {
                o(e), h(!0);
              },
              refresh: j,
            }),
            Object(xe.jsx)(s.LoadingIndicatorDialog, {
              handleClose: function () {
                v();
              },
              area: "global_loading_indicator",
            }),
            Object(xe.jsx)(ze, {
              openSDFDialog: d,
              handleClose: function (e) {
                h(!1),
                  null !== e &&
                    ((x = new AbortController()),
                    new Me().resolvePath(a, t, O, e, x));
              },
            }),
          ],
        });
      }
      var Xe = n(273),
        Ze = n(258);
      function $e(e) {
        var t = e.children;
        return Object(xe.jsx)(Ze.a, {
          variant: "outlined",
          position: "relative",
          color: "transparent",
          children: Object(xe.jsx)(Xe.a, { children: t }),
        });
      }
      function et() {
        return Object(xe.jsxs)($e, {
          children: [
            Object(xe.jsx)("a", {
              href: "https://www.bayer.com",
              target: "_blank",
              children: Object(xe.jsx)("img", {
                style: { height: 48, marginLeft: 48 },
                src: "textures/bayer-logo.svg",
                alt: "Powered By Bayer",
              }),
            }),
            Object(xe.jsx)(S.a, {
              variant: "h6",
              style: { marginLeft: 48, color: "rgba(0, 0, 0, 0.54)" },
              children: "CIME: Chem-Informatics Model Explorer",
            }),
          ],
        });
      }
      var tt = !1;
      s.PluginRegistry.getInstance().registerPlugin(new Te());
      var nt = function () {
        var e = Object(r.useState)(new s.API()),
          t = Object(c.a)(e, 2),
          n = t[0];
        return (
          t[1],
          Object(xe.jsx)(s.PSEContextProvider, {
            context: n,
            children: Object(xe.jsx)(s.Application, {
              config: { preselect: { url: "datasets/test.sdf" } },
              features: { disableEmbeddings: { tsne: !0, forceatlas: !0 } },
              overrideComponents: { datasetTab: Qe, appBar: et },
            }),
          })
        );
      };
      o.a.render(
        Object(xe.jsx)(i.a.StrictMode, { children: Object(xe.jsx)(nt, {}) }),
        document.getElementById("root")
      );
    },
  },
  [[342, 1, 2]],
]);
//# sourceMappingURL=main.8a28a57a.chunk.js.map