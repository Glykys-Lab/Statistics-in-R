function saveWorksheetCell(event) {
  const cell = event.target.closest('.worksheet-cell');
  if (!cell || !window.Shiny) return;
  Shiny.setInputValue('cell_edit', {
    row: Number(cell.dataset.row), column: Number(cell.dataset.column),
    value: cell.value, revision: Number(cell.dataset.revision), nonce: Date.now()
  }, {priority: 'event'});
}
// Capture before Shiny's ordinary text-input bindings handle the same event.
document.addEventListener('change', saveWorksheetCell, true);
document.addEventListener('blur', saveWorksheetCell, true);
document.addEventListener('click', function(event) {
  const link = event.target.closest('.project-link');
  if (link && window.Shiny) Shiny.setInputValue('project_open', {
    kind: link.dataset.kind, key: link.dataset.key, nonce: Date.now()
  }, {priority: 'event'});
}, true);
// Keyboard movement within the worksheet keeps the familiar spreadsheet order.
document.addEventListener('keydown', function (event) {
  const cell = event.target.closest('.worksheet-cell');
  if (!cell || !['ArrowUp', 'ArrowDown', 'Enter'].includes(event.key)) return;
  event.preventDefault();
  const row = Number(cell.dataset.row) + (event.key === 'ArrowUp' ? -1 : 1);
  const next = document.querySelector('.worksheet-cell[data-row="' + row + '"][data-column="' + cell.dataset.column + '"]');
  if (next) { next.focus(); next.select(); }
});
