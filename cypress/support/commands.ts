/// <reference types="cypress" />
// ***********************************************
// This example commands.ts shows you how to
// create various custom commands and overwrite
// existing commands.
//
// For more comprehensive examples of custom
// commands please read more here:
// https://on.cypress.io/custom-commands
// ***********************************************
//
//
// -- This is a parent command --
// Cypress.Commands.add('login', (email, password) => { ... })
//
//
// -- This is a child command --
// Cypress.Commands.add('drag', { prevSubject: 'element'}, (subject, options) => { ... })
//
//
// -- This is a dual command --
// Cypress.Commands.add('dismiss', { prevSubject: 'optional'}, (subject, options) => { ... })
//
//
// -- This will overwrite an existing command --
// Cypress.Commands.overwrite('visit', (originalFn, url, options) => { ... })
//
// declare global {
//   namespace Cypress {
//     interface Chainable {
//       login(email: string, password: string): Chainable<void>
//       drag(subject: string, options?: Partial<TypeOptions>): Chainable<Element>
//       dismiss(subject: string, options?: Partial<TypeOptions>): Chainable<Element>
//       visit(originalFn: CommandOriginalFn, url: string, options: Partial<VisitOptions>): Chainable<Element>
//     }
//   }
// }

Cypress.Commands.add('createProject', (name: string, path: string) => {
  // Intercept the network request
  cy.intercept('POST', '/api/reaction_cime/upload_csv').as('loginMutation');

  cy.byId('upload-file-input').selectFile(path, { force: true }).wait(500);

  // eslint-disable-next-line cypress/no-unnecessary-waiting
  cy.wait('@loginMutation');
});

Cypress.Commands.add('selectSomething', () => {
  // eslint-disable-next-line cypress/no-unnecessary-waiting
  cy.get('[data-cy="small-multiple-container"]')
    .wait(1000)
    .trigger('mousedown', 100, 100, { button: 0, offsetX: 100, offsetY: 100 })
    .wait(1000)
    .trigger('mousemove', 600, 100, { button: 0, offsetX: 600, offsetY: 100 })
    .wait(1000)
    .trigger('mousemove', 600, 600, { button: 0, offsetX: 600, offsetY: 600 })
    .wait(1000)
    .trigger('mousemove', 100, 600, { button: 0, offsetX: 100, offsetY: 600 })
    .wait(1000)
    .trigger('mouseup', 100, 600, { button: 0, offsetX: 100, offsetY: 600 })
    .wait(1000)
    .rightclick();
});

Cypress.Commands.add('deleteProject', () => {
  cy.get('[data-cy="delete-project-card-action"]').click();
  cy.get('[data-cy="suredelete"]').click();
});

Cypress.Commands.add('saveSessionAs', (name: string) => {
  cy.get('[data-cy="save-session-as-context-button"]').click();
  // eslint-disable-next-line cypress/no-unnecessary-waiting
  cy.get('[data-cy="cimeSaveSessionName"]').wait(1000).focus().type(name);
  // eslint-disable-next-line cypress/no-unnecessary-waiting
  cy.get('[data-cy="savesessionbtn"]').click();
});

Cypress.Commands.add('switchTab', (name: 'dataset' | 'encoding' | 'projection') => {
  cy.get(`[data-cy="${name}-tab"]`).click();
});

Cypress.Commands.add('chooseNth', (nth: number) => {
  cy.get(`.MuiAutocomplete-popper li[data-option-index="${nth}"]`).click();
});

Cypress.Commands.add('byId', (dataCy: string) => {
  return cy.get(`[data-cy=${dataCy}]`);
});

Cypress.Commands.add('URLParams', (url: string) => {
  const arr = url.split('/?')[1];

  const paramObj = {};
  if (arr) {
    arr.split('&').forEach((param) => {
      const [key, value] = param.split('=');
      paramObj[key] = value;
    });
  }

  return cy.wrap(paramObj);
});
