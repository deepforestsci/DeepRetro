## UI debugging

- Reproduce responsive layout bugs with the user's exact item count and configuration. A one-item fixture does not validate multi-item sidebar sizing, scroll containment, or grid overflow.
- Before changing responsive CSS, inspect `scrollWidth`, `clientWidth`, and the bounding boxes of the shell, workspace, columns, and overflowing descendants in the browser.
